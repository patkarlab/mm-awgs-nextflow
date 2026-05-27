#!/usr/bin/env python3
"""
build_v6_panel.py
=================

Build the v6 MM-focused adaptive sampling panel from v5.

What this does:
  1. Reads the v5 panel BED (T2T, chr-naming) verbatim.
  2. Removes the v5 TP53 line (will be rebuilt with a tighter flank).
  3. Looks up coordinates for TP53 + 10 new MM driver genes from the
     supplied T2T RefSeq GFF (NC_-named contigs).
  4. Applies per-gene flanks, clips to contig lengths from the T2T .fai.
  5. Merges any newly-introduced overlaps, sorts in natural chr order.
  6. Writes the two BEDs that downstream tools actually consume:
       aWGS_MMfocused_v6_t2t_chr.bed   (for analysis on gandalf)
       aWGS_MMfocused_v6_t2t_NC.bed    (for MinKNOW on the P2i)
     Both are the same coordinates; only the contig naming differs.

Design notes
------------
- Variant-agnostic and finding-agnostic. The only hardcoded data is:
    - 10 gene symbols (DIS3, TRAF3, PRDM1, ATM, CYLD, H1-4, MAX, EGR1, LTB, ATR)
    - per-gene flanks (all +/-50 kb except TP53 which is +/-500 kb)
    - TP53's symbol (for the replace-in-place trim)
    - HGNC alias map (H1-4 / HIST1H1E only)
  All coordinates come from the annotation file at runtime.
- v5 retained regions are NEVER re-derived. They are copied byte-for-byte
  from the v5 BED supplied at runtime, so any region the lab already
  validated keeps its exact coordinates in v6.
- TP53 is identified for replacement by its BED label exactly ("TP53").

Usage
-----
  python3 build_v6_panel.py \\
    --v5-bed   /goast/nikhil_awgs_testing/t2t/beds/aWGS_MMfocused_v5_t2t_chr.bed \\
    --t2t-gff  /goast/hemat_data/references/T2T/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff \\
    --t2t-fai  /goast/nikhil_awgs_testing/t2t/refs/chm13v2.0.ucsc.fa.fai \\
    --outdir   /goast/nikhil_awgs_testing/panel/v6_build/

The script writes nothing outside --outdir.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import logging
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ---------------------------------------------------------------------------
# v6 design (genes and flanks only; coordinates come from annotation)
# ---------------------------------------------------------------------------

# Label used in the v5 BED for the region we are replacing.
TP53_LABEL = "TP53"
TP53_FLANK_BP = 500_000   # v6 trim: was +/-1 Mb in v5

# 10 new MM driver genes for v6, each as a focal region (+/-50 kb).
# primary_symbol -> (flank_bp, list_of_aliases_tried_in_order_if_primary_misses)
NEW_GENES: Dict[str, Tuple[int, List[str]]] = {
    "DIS3":  (50_000, []),
    "TRAF3": (50_000, []),
    "PRDM1": (50_000, []),
    "ATM":   (50_000, []),
    "CYLD":  (50_000, []),
    "H1-4":  (50_000, ["HIST1H1E"]),  # current HGNC + legacy alias
    "MAX":   (50_000, []),
    "EGR1":  (50_000, []),
    "LTB":   (50_000, []),
    "ATR":   (50_000, []),
}

# T2T-CHM13v2.0 NCBI RefSeq accession <-> chr-name mapping (autosomes + sex).
NC_TO_CHR: Dict[str, str] = {
    "NC_060925.1": "chr1",  "NC_060926.1": "chr2",  "NC_060927.1": "chr3",
    "NC_060928.1": "chr4",  "NC_060929.1": "chr5",  "NC_060930.1": "chr6",
    "NC_060931.1": "chr7",  "NC_060932.1": "chr8",  "NC_060933.1": "chr9",
    "NC_060934.1": "chr10", "NC_060935.1": "chr11", "NC_060936.1": "chr12",
    "NC_060937.1": "chr13", "NC_060938.1": "chr14", "NC_060939.1": "chr15",
    "NC_060940.1": "chr16", "NC_060941.1": "chr17", "NC_060942.1": "chr18",
    "NC_060943.1": "chr19", "NC_060944.1": "chr20", "NC_060945.1": "chr21",
    "NC_060946.1": "chr22", "NC_060947.1": "chrX",  "NC_060948.1": "chrY",
}
CHR_TO_NC: Dict[str, str] = {v: k for k, v in NC_TO_CHR.items()}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class Region:
    """A single panel region in chr-style coordinates (BED 0-based half-open)."""
    chrom: str
    start: int
    end: int
    name: str
    provenance: str   # 'v5_retained' | 'v5_modified' | 'v6_new' | 'merged'
    flank_bp: Optional[int] = None
    raw_gene_start: Optional[int] = None
    raw_gene_end: Optional[int] = None


@dataclass
class GeneCoord:
    """Whole-gene span as discovered from the annotation."""
    chrom: str   # NC_-named (from the T2T GFF)
    start: int   # 0-based, half-open
    end: int     # exclusive
    symbol: str  # symbol as written in the annotation


# ---------------------------------------------------------------------------
# GFF parsing
# ---------------------------------------------------------------------------

def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def _parse_gff3_attrs(field: str) -> Dict[str, str]:
    """GFF3 column 9: key=value;key=value;..."""
    out: Dict[str, str] = {}
    for chunk in field.strip().split(";"):
        chunk = chunk.strip()
        if not chunk or "=" not in chunk:
            continue
        k, v = chunk.split("=", 1)
        out[k] = v
    return out


def parse_gff(path: Path, wanted: set) -> Dict[str, GeneCoord]:
    """
    Whole-gene span lookup from an NCBI RefSeq GFF3.

    Returns a dict keyed by uppercase gene symbol. If a symbol appears on
    more than one contig (e.g. an alt placement), the NC_ primary contig
    is preferred.
    """
    wanted_u = {s.upper() for s in wanted}
    span = defaultdict(lambda: {"chrom": None, "start": None, "end": None, "symbol": None})
    with _open_text(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = _parse_gff3_attrs(cols[8])
            sym = attrs.get("Name") or attrs.get("gene") or ""
            sym_u = sym.upper()
            if sym_u not in wanted_u:
                continue
            chrom = cols[0]
            start0 = int(cols[3]) - 1
            end1 = int(cols[4])
            rec = span[sym_u]
            if rec["chrom"] is None:
                rec["chrom"] = chrom
                rec["symbol"] = sym
            elif rec["chrom"] != chrom:
                # Prefer NC_ primary contig over any alt/scaffold.
                if chrom.startswith("NC_") and not rec["chrom"].startswith("NC_"):
                    rec["chrom"] = chrom
                    rec["start"] = None
                    rec["end"] = None
                else:
                    continue
            if rec["start"] is None or start0 < rec["start"]:
                rec["start"] = start0
            if rec["end"] is None or end1 > rec["end"]:
                rec["end"] = end1
    return {
        u: GeneCoord(r["chrom"], r["start"], r["end"], r["symbol"])
        for u, r in span.items()
        if None not in (r["chrom"], r["start"], r["end"])
    }


def lookup(coords: Dict[str, GeneCoord], primary: str, aliases: List[str]) -> Optional[GeneCoord]:
    for sym in [primary] + aliases:
        if sym.upper() in coords:
            return coords[sym.upper()]
    return None


def report_missing(coords: Dict[str, GeneCoord], log: logging.Logger) -> List[str]:
    """List primary symbols (incl. TP53) that had no annotation hit."""
    missing: List[str] = []
    if lookup(coords, "TP53", []) is None:
        missing.append("TP53")
    for primary, (_f, aliases) in NEW_GENES.items():
        if lookup(coords, primary, aliases) is None:
            missing.append(primary)
    if missing:
        log.error("Annotation lookup failed for: %s", missing)
    return missing


# ---------------------------------------------------------------------------
# FAI and BED I/O
# ---------------------------------------------------------------------------

def parse_fai(path: Path) -> Dict[str, int]:
    """Read a samtools .fai index into {contig_name: length_bp}."""
    sizes: Dict[str, int] = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes


def read_bed_as_regions(path: Path) -> List[Region]:
    """Read a 4-column BED into Region objects, provenance 'v5_retained'."""
    regs: List[Region] = []
    with open(path) as fh:
        for ln, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                logging.warning("%s line %d: skipping malformed BED row", path, ln)
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            regs.append(Region(chrom=chrom, start=start, end=end,
                               name=name, provenance="v5_retained"))
    return regs


def write_bed(
    regions: List[Region],
    path: Path,
    chrom_map: Optional[Dict[str, str]] = None,
) -> int:
    """
    Write a 4-column BED, sorted by natural chrom order then start.

    If chrom_map is supplied, rows whose chrom is not a key in chrom_map
    are dropped with a warning (so an inconsistent BED is never written).
    """
    sorted_regs = sorted(regions, key=lambda r: (_chrom_sort_key(r.chrom), r.start))
    n_written = 0
    missing: List[str] = []
    with open(path, "w") as fh:
        for r in sorted_regs:
            if chrom_map is not None:
                if r.chrom not in chrom_map:
                    missing.append(r.chrom)
                    continue
                out_chrom = chrom_map[r.chrom]
            else:
                out_chrom = r.chrom
            fh.write(f"{out_chrom}\t{r.start}\t{r.end}\t{r.name}\n")
            n_written += 1
    if missing:
        logging.error("write_bed(%s): %d rows dropped (chrom not in map): %s",
                      path.name, len(missing), sorted(set(missing)))
    return n_written


# ---------------------------------------------------------------------------
# Region operations
# ---------------------------------------------------------------------------

def apply_flank_clip(
    chrom: str,
    gene_start: int,
    gene_end: int,
    flank: int,
    sizes: Dict[str, int],
) -> Tuple[int, int]:
    """Apply a symmetric flank and clip to [0, contig_length]."""
    start = max(0, gene_start - flank)
    end = gene_end + flank
    if chrom in sizes:
        end = min(end, sizes[chrom])
    return start, end


def _chrom_sort_key(c: str) -> Tuple[int, int, str]:
    """Natural-numeric ordering: chr1..chr22, chrX, chrY, others alpha."""
    stripped = c[3:] if c.startswith("chr") else c
    if stripped.isdigit():
        return (0, int(stripped), "")
    if stripped == "X":
        return (1, 0, "")
    if stripped == "Y":
        return (2, 0, "")
    return (3, 0, stripped)


def merge_overlaps(regions: List[Region]) -> List[Region]:
    """Sort by (chrom, start) and merge overlapping or touching regions."""
    if not regions:
        return []
    regs = sorted(regions, key=lambda r: (_chrom_sort_key(r.chrom), r.start))
    out: List[Region] = [regs[0]]
    for r in regs[1:]:
        last = out[-1]
        if r.chrom == last.chrom and r.start <= last.end:
            existing = last.name.split("+")
            new_name = last.name if r.name in existing else f"{last.name}+{r.name}"
            out[-1] = Region(
                chrom=last.chrom,
                start=last.start,
                end=max(last.end, r.end),
                name=new_name,
                provenance="merged",
            )
        else:
            out.append(r)
    return out


def drop_invalid(regions: List[Region]) -> Tuple[List[Region], List[Region]]:
    """Split into (valid, invalid). Invalid means end <= start."""
    valid: List[Region] = []
    invalid: List[Region] = []
    for r in regions:
        if r.end <= r.start:
            logging.error("Dropping zero/negative-width region %s:%d-%d (%s)",
                          r.chrom, r.start, r.end, r.name)
            invalid.append(r)
        else:
            valid.append(r)
    return valid, invalid


def detect_overlaps(regions: List[Region]) -> List[str]:
    """After merge, no overlaps should exist; report any that slipped through."""
    issues: List[str] = []
    prev: Optional[Region] = None
    for r in sorted(regions, key=lambda x: (_chrom_sort_key(x.chrom), x.start)):
        if prev is not None and prev.chrom == r.chrom and r.start < prev.end:
            issues.append(
                f"residual overlap {prev.chrom}:{prev.start}-{prev.end} ({prev.name}) "
                f"vs {r.chrom}:{r.start}-{r.end} ({r.name})"
            )
        prev = r
    return issues


def _sizes_chr_keyed(sizes: Dict[str, int]) -> Dict[str, int]:
    """Return a copy of a sizes dict translating any NC_ keys to chr names."""
    out: Dict[str, int] = {}
    for k, v in sizes.items():
        if k.startswith("NC_") and k in NC_TO_CHR:
            out[NC_TO_CHR[k]] = v
        else:
            out[k] = v
    return out


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def build_panel(
    v5_regions: List[Region],
    coords: Dict[str, GeneCoord],
    sizes: Dict[str, int],
) -> Tuple[List[Region], List[Region], List[str]]:
    """
    Build v6 regions from v5 + annotation.

    Returns: (final_regions_sorted_and_merged, invalid_dropped, residual_overlap_warnings)
    Coordinates in 'coords' are NC_-named; this function translates to chr-names.
    """
    # 1. Start from v5, removing the TP53 line.
    retained = [r for r in v5_regions if r.name != TP53_LABEL]
    n_removed = len(v5_regions) - len(retained)
    logging.info("Removed %d v5 line(s) labeled '%s'; keeping %d v5 regions",
                 n_removed, TP53_LABEL, len(retained))
    if n_removed != 1:
        logging.warning("Expected exactly 1 TP53 line in v5; found %d", n_removed)

    new_regs: List[Region] = []
    sizes_chr = _sizes_chr_keyed(sizes)

    # 2a. TP53 (modified flank).
    tp53 = lookup(coords, "TP53", [])
    if tp53 is not None:
        chrom = NC_TO_CHR.get(tp53.chrom, tp53.chrom)
        s, e = apply_flank_clip(chrom, tp53.start, tp53.end, TP53_FLANK_BP, sizes_chr)
        new_regs.append(Region(
            chrom=chrom, start=s, end=e, name=TP53_LABEL,
            provenance="v5_modified", flank_bp=TP53_FLANK_BP,
            raw_gene_start=tp53.start, raw_gene_end=tp53.end,
        ))

    # 2b. 10 new genes.
    for primary, (flank, aliases) in NEW_GENES.items():
        gc = lookup(coords, primary, aliases)
        if gc is None:
            continue  # error already logged by report_missing()
        chrom = NC_TO_CHR.get(gc.chrom, gc.chrom)
        s, e = apply_flank_clip(chrom, gc.start, gc.end, flank, sizes_chr)
        new_regs.append(Region(
            chrom=chrom, start=s, end=e, name=primary,
            provenance="v6_new", flank_bp=flank,
            raw_gene_start=gc.start, raw_gene_end=gc.end,
        ))

    # 3. Combine, drop invalid, merge overlaps.
    combined = retained + new_regs
    valid, invalid = drop_invalid(combined)
    merged = merge_overlaps(valid)
    residual = detect_overlaps(merged)
    return merged, invalid, residual


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------

def total_bp(regs: List[Region]) -> int:
    return sum(r.end - r.start for r in regs)


def md5_of(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def write_build_report(regions: List[Region], out_path: Path) -> None:
    """TSV documenting every region in the v6 panel."""
    with open(out_path, "w") as fh:
        fh.write("\t".join([
            "name", "chrom", "start", "end", "size_bp",
            "provenance", "flank_bp",
            "raw_gene_start", "raw_gene_end",
        ]) + "\n")
        for r in sorted(regions, key=lambda x: (_chrom_sort_key(x.chrom), x.start)):
            fh.write("\t".join([
                r.name, r.chrom, str(r.start), str(r.end), str(r.end - r.start),
                r.provenance,
                str(r.flank_bp) if r.flank_bp is not None else "NA",
                str(r.raw_gene_start) if r.raw_gene_start is not None else "NA",
                str(r.raw_gene_end) if r.raw_gene_end is not None else "NA",
            ]) + "\n")


def write_stamp(
    outdir: Path,
    bed_paths: List[Path],
    inputs: Dict[str, Path],
    invalid: List[Region],
    residual: List[str],
    missing: List[str],
) -> None:
    stamp = outdir / "aWGS_MMfocused_v6_BUILD_STAMP.txt"
    with open(stamp, "w") as fh:
        fh.write("Panel v6 build stamp\n")
        fh.write("====================\n\n")
        fh.write("Inputs:\n")
        for k, v in inputs.items():
            fh.write(f"  {k}: {v}\n")
        fh.write("\nOutputs:\n")
        for bed in bed_paths:
            with open(bed) as bfh:
                lines = [l for l in bfh if l.strip()]
            n_regions = len(lines)
            bp = sum(int(l.split('\t')[2]) - int(l.split('\t')[1]) for l in lines)
            fh.write(f"  {bed.name}\n")
            fh.write(f"    MD5={md5_of(bed)}\n")
            fh.write(f"    regions={n_regions}\n")
            fh.write(f"    size_bp={bp}\n")
            fh.write(f"    size_Mb={bp / 1e6:.2f}\n")
        fh.write("\nGenes missing from annotation:\n")
        fh.write("  (none)\n" if not missing else "".join(f"  {g}\n" for g in missing))
        fh.write("\nDropped invalid regions (end <= start):\n")
        fh.write("  (none)\n" if not invalid else
                 "".join(f"  {r.chrom}:{r.start}-{r.end}  {r.name}\n" for r in invalid))
        fh.write("\nResidual overlap warnings:\n")
        fh.write("  (none)\n" if not residual else
                 "".join(f"  {s}\n" for s in residual))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--v5-bed", required=True, type=Path,
                    help="v5 panel BED in chr-style T2T coords (4-column).")
    ap.add_argument("--t2t-gff", required=True, type=Path,
                    help="NCBI RefSeq T2T-CHM13v2.0 GFF (NC_-named).")
    ap.add_argument("--t2t-fai", required=True, type=Path,
                    help=".fai for the chr-named T2T FASTA used in analysis.")
    ap.add_argument("--outdir", required=True, type=Path,
                    help="Output directory (created if needed).")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    log = logging.getLogger("build_v6")
    args.outdir.mkdir(parents=True, exist_ok=True)

    # Build the lookup set: TP53 + 10 new genes + their aliases.
    wanted: set = {"TP53"}
    for primary, (_f, aliases) in NEW_GENES.items():
        wanted.add(primary)
        wanted.update(aliases)
    log.info("Looking up %d gene symbols (incl. aliases): %s",
             len(wanted), sorted(wanted))

    log.info("Parsing T2T GFF: %s", args.t2t_gff)
    coords = parse_gff(args.t2t_gff, wanted)
    missing = report_missing(coords, log)

    log.info("Reading v5 BED: %s", args.v5_bed)
    v5_regions = read_bed_as_regions(args.v5_bed)
    log.info("v5: %d regions, %.2f Mb", len(v5_regions), total_bp(v5_regions) / 1e6)

    log.info("Parsing T2T FAI: %s", args.t2t_fai)
    t2t_sizes = parse_fai(args.t2t_fai)

    log.info("Building v6 panel")
    final_regs, invalid, residual = build_panel(v5_regions, coords, t2t_sizes)

    # Write the two BEDs.
    bed_chr = args.outdir / "aWGS_MMfocused_v6_t2t_chr.bed"
    bed_nc  = args.outdir / "aWGS_MMfocused_v6_t2t_NC.bed"
    n1 = write_bed(final_regs, bed_chr)
    n2 = write_bed(final_regs, bed_nc, chrom_map=CHR_TO_NC)

    log.info("v6 chr BED: %d regions, %.2f Mb", n1, total_bp(final_regs) / 1e6)
    log.info("v6 NC  BED: %d regions", n2)

    if n1 != n2:
        log.error("chr and NC BEDs have different region counts (%d vs %d). "
                  "A chrom is missing from NC_TO_CHR.", n1, n2)

    # Build report and stamp.
    report = args.outdir / "aWGS_MMfocused_v6_build_report.tsv"
    write_build_report(final_regs, report)
    write_stamp(
        args.outdir,
        [bed_chr, bed_nc],
        {
            "v5 BED":  args.v5_bed,
            "T2T GFF": args.t2t_gff,
            "T2T FAI": args.t2t_fai,
        },
        invalid, residual, missing,
    )
    log.info("Wrote: %s", report)
    log.info("Wrote: %s", args.outdir / "aWGS_MMfocused_v6_BUILD_STAMP.txt")

    if missing or invalid or residual or n1 != n2:
        log.error("Build completed with warnings; review the stamp file before use.")
        return 2
    log.info("Build OK.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
