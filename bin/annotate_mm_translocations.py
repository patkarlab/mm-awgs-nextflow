#!/usr/bin/env python3
"""
annotate_mm_translocations.py
=============================

Annotates a SURVIVOR-merged VCF of structural variants against the MM-specific
translocation panel by intersecting BND breakpoints with the panel BED.

For each record where at least one breakpoint falls inside a panel region,
emits a row of:
  - sample id
  - SV id / type / filter
  - both breakpoints (chrom, pos, gene region)
  - whether the pair matches a known MM partner pair (from dictionary)
  - supporting callers inferred from SURVIVOR's SUPP_VEC
  - support read counts from FORMAT/DV (fallback INFO RE / SUPPORT / SR)

The dictionary file is the only source of biological priors. The script
itself never hardcodes any expected breakpoints, sample-specific findings,
or known karyotypes.

Usage:
  annotate_mm_translocations.py \\
      --vcf merged.vcf.gz \\
      --panel-bed aWGS_MMfocused_v6_t2t_chr.bed \\
      --dictionary mm_translocation_dictionary.tsv \\
      --sample SAMPLE_ID \\
      --output SAMPLE_ID.mm_annotated.tsv
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


__version__ = "0.3.0"
# [dictionary-token-matching applied]
# [cytoband-partner-annotation applied]


@dataclass
class PanelRegion:
    chrom: str
    start: int
    end: int
    name: str

    def contains(self, chrom: str, pos: int) -> bool:
        return chrom == self.chrom and self.start <= pos < self.end


@dataclass
class SvRecord:
    chrom: str
    pos: int
    sv_id: str
    sv_type: str
    mate_chrom: Optional[str]
    mate_pos: Optional[int]
    filt: str
    info: Dict[str, str]
    callers: List[str]
    support: str = ""


def load_panel(bed_path: Path) -> List[PanelRegion]:
    out: List[PanelRegion] = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            out.append(PanelRegion(chrom, start, end, name))
    return out


# Region-suffix words that are not gene symbols (e.g. the "locus" in
# "IGH_locus"). Dropped during tokenization so IGH_locus and IGH match.
_PARTNER_SUFFIX_TOKENS = {"LOCUS"}


def _norm_tokens(name: str) -> frozenset:
    """Reduce a partner label to a set of normalized gene tokens.

    Splits on '/' and '_' (and '+'), uppercases, and drops region-suffix
    words. E.g. 'IGH_locus' -> {IGH}; 'FGFR3/NSD2' -> {FGFR3, NSD2};
    'IGL/IGLL5' -> {IGL, IGLL5}. Two partners are considered the same locus
    when their token sets intersect.
    """
    raw = str(name).strip().upper().replace("+", "/").replace("_", "/")
    return frozenset(t for t in raw.split("/")
                     if t and t not in _PARTNER_SUFFIX_TOKENS)


def load_dictionary(dict_path):
    """
    Load the MM translocation dictionary as a list of
    (tokset_a, tokset_b, row) entries. Partner names are reduced to gene
    token sets (see _norm_tokens) so lookup matches on gene identity rather
    than exact label. Missing dictionary file is non-fatal (empty list).
    """
    out = []
    if not dict_path.exists():
        return out
    with open(dict_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            a = (row.get("partner_a") or "").strip()
            b = (row.get("partner_b") or "").strip()
            if not a or not b:
                continue
            out.append((_norm_tokens(a), _norm_tokens(b), row))
    return out


def dictionary_lookup(dictionary, gene_a, gene_b):
    """Return the dictionary row for an unordered gene-a / gene-b pair, or
    None. Matches when the pair's token sets intersect the entry's token
    sets on both sides (in either orientation)."""
    ta, tb = _norm_tokens(gene_a), _norm_tokens(gene_b)
    if not ta or not tb:
        return None
    for da, db, row in dictionary:
        if (ta & da and tb & db) or (ta & db and tb & da):
            return row
    return None


def parse_info(info_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in info_field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = "True"
    return out


def parse_bnd_alt(alt: str) -> Tuple[Optional[str], Optional[int]]:
    """Parse BND ALT (N]chr14:99000000] or [chr8:130000000[N)."""
    for opener, closer in [("]", "]"), ("[", "[")]:
        if opener in alt:
            try:
                inner = alt.split(opener)[1].split(closer)[0]
                chrom, pos = inner.rsplit(":", 1)
                return chrom, int(pos)
            except (IndexError, ValueError):
                continue
    return None, None


def infer_callers(info: Dict[str, str]) -> List[str]:
    """
    Decode SUPP_VEC bitstring written by SURVIVOR. Bit order matches the
    order of input VCFs passed to SURVIVOR merge; in our pipeline that is
    [Sniffles, CuteSV, Severus].
    """
    out = []
    supp_vec = info.get("SUPP_VEC", "")
    caller_order = ["Sniffles", "CuteSV", "Severus"]
    for bit, name in zip(supp_vec, caller_order):
        if bit == "1":
            out.append(name)
    return out


def support_reads_from(info: Dict[str, str], fmt: str, sample_cols: List[str]) -> str:
    """
    Variant-supporting read count for a record, read generically (no hardcoded
    calls). The merged VCF carries one sample column per input caller; the
    caller that detected the junction holds the real FORMAT/DV while the others
    are zero/absent, so DV is taken as the MAX across all sample columns. If no
    DV is present, fall back to caller-specific INFO tags:
      RE       (CuteSV)
      SUPPORT  (Sniffles)
      SR       (generic)
      SUPP_READS first colon-field (Severus, e.g. "2:0:0:2:0:0" -> 2)
    """
    if fmt:
        keys = fmt.split(":")
        best = None
        for col in sample_cols:
            if not col or col in (".", "./."):
                continue
            f = dict(zip(keys, col.split(":")))
            dv = f.get("DV", "")
            if dv.isdigit():
                best = int(dv) if best is None else max(best, int(dv))
        if best is not None:
            return str(best)
    for tag in ("RE", "SUPPORT", "SR"):
        v = str(info.get(tag, ""))
        if v.isdigit():
            return v
    sr = str(info.get("SUPP_READS", ""))
    first = sr.split(":")[0] if sr else ""
    if first.isdigit():
        return first
    return ""


def open_vcf(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def parse_vcf(vcf_path: Path) -> List[SvRecord]:
    out: List[SvRecord] = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, sv_id, _ref, alt, _qual, filt, info_field = cols[:8]
            pos = int(pos)
            info = parse_info(info_field)
            sv_type = info.get("SVTYPE", "")
            mate_chrom, mate_pos = (None, None)

            if sv_type in ("BND", "TRA"):
                # BND records have mate coords in ALT (e.g. N]chr11:69500000]).
                # SURVIVOR renames cross-chromosome BNDs to TRA while keeping
                # the same ALT format, so we parse both the same way.
                mate_chrom, mate_pos = parse_bnd_alt(alt)
                # Fallback: if ALT parsing failed (some SURVIVOR records use
                # symbolic ALT like <TRA>), use SURVIVOR's CHR2/END INFO tags.
                if mate_chrom is None:
                    chr2 = info.get("CHR2")
                    end = info.get("END")
                    if chr2 and end:
                        try:
                            mate_chrom, mate_pos = chr2, int(end)
                        except ValueError:
                            pass

            elif sv_type in ("DEL", "DUP", "INV", "INS"):
                # Same-chromosome events: mate side is END on the same chrom.
                end = info.get("END")
                if end:
                    try:
                        mate_chrom, mate_pos = chrom, int(end)
                    except ValueError:
                        pass

            # FORMAT + per-caller sample columns carry read support (DV).
            fmt = cols[8] if len(cols) > 8 else ""
            sample_cols = cols[9:] if len(cols) > 9 else []
            support = support_reads_from(info, fmt, sample_cols)

            out.append(SvRecord(
                chrom=chrom, pos=pos, sv_id=sv_id, sv_type=sv_type,
                mate_chrom=mate_chrom, mate_pos=mate_pos,
                filt=filt, info=info, callers=infer_callers(info),
                support=support,
            ))
    return out



@dataclass
class CytobandTable:
    """Cytoband lookup for T2T-CHM13v2.0. Bands tile each chromosome
    contiguously, so any in-range coordinate resolves to exactly one band."""
    # {chrom: [(start, end, band), ...] sorted by start}
    bands: Dict[str, List[Tuple[int, int, str]]]

    def band_for(self, chrom: Optional[str], pos: Optional[int]) -> Optional[str]:
        if chrom is None or pos is None:
            return None
        arr = self.bands.get(chrom)
        if not arr:
            return None
        # Linear scan is fine: <= ~60 bands per chromosome.
        for start, end, name in arr:
            if start <= pos < end:
                return name
        return None


def load_cytobands(bed_path: Path) -> CytobandTable:
    """Load a 5-column UCSC cytoband BED (chrom start end band gieStain).
    Only the first four columns are used."""
    bands: Dict[str, List[Tuple[int, int, str]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, band = parts[0], int(parts[1]), int(parts[2]), parts[3]
            bands.setdefault(chrom, []).append((start, end, band))
    for chrom in bands:
        bands[chrom].sort(key=lambda t: t[0])
    if not bands:
        sys.stderr.write(f"ERROR: no cytobands parsed from {bed_path}\n")
        sys.exit(1)
    return CytobandTable(bands)


def _strip_chr(chrom: str) -> str:
    """chr8 -> 8, chrX -> X, for clinical band notation (8q24.21)."""
    return chrom[3:] if chrom.startswith("chr") else chrom


def characterize_side(chrom, pos, region, cytobands):
    """Return (label, source) for one breakpoint side.

    panel region      -> (region.name, "panel")
    else band found   -> ("8q24.21",   "cytoband")
    else (defensive)  -> ("chr8:127.0Mb", "coordinate")
    """
    if region is not None:
        return region.name, "panel"
    band = cytobands.band_for(chrom, pos)
    if band is not None:
        return f"{_strip_chr(chrom)}{band}", "cytoband"
    if chrom is not None and pos is not None:
        return f"{chrom}:{pos / 1e6:.1f}Mb", "coordinate"
    return "OFF_PANEL", "coordinate"


def region_for(chrom: Optional[str], pos: Optional[int], panel: List[PanelRegion]) -> Optional[PanelRegion]:
    if chrom is None or pos is None:
        return None
    for r in panel:
        if r.contains(chrom, pos):
            return r
    return None


def annotate(records, panel, dictionary, sample, cytobands):
    out = []
    for r in records:
        side_a = region_for(r.chrom, r.pos, panel)
        side_b = region_for(r.mate_chrom, r.mate_pos, panel)
        if side_a is None and side_b is None:
            continue

        gene_a, gene_a_source = characterize_side(r.chrom, r.pos, side_a, cytobands)
        gene_b, gene_b_source = characterize_side(r.mate_chrom, r.mate_pos, side_b, cytobands)

        known = ""
        freq = ""
        if side_a and side_b and gene_a_source == "panel" and gene_b_source == "panel":
            hit = dictionary_lookup(dictionary, gene_a, gene_b)
            if hit:
                known = hit.get("name", "yes")
                freq = hit.get("frequency", "")

        out.append({
            "sample":         sample,
            "sv_id":          r.sv_id,
            "sv_type":        r.sv_type,
            "filter":         r.filt,
            "chrom_a":        r.chrom,
            "pos_a":          str(r.pos),
            "gene_a":         gene_a,
            "chrom_b":        r.mate_chrom or "",
            "pos_b":          str(r.mate_pos) if r.mate_pos is not None else "",
            "gene_b":         gene_b,
            "gene_a_source":  gene_a_source,
            "gene_b_source":  gene_b_source,
            "known_mm_pair":  known,
            "known_freq":     freq,
            "callers":        ",".join(r.callers) or "unknown",
            "n_callers":      str(len(r.callers)),
            "supp_vec":       r.info.get("SUPP_VEC", ""),
            "support_reads":  r.support,
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf",         required=True, type=Path)
    ap.add_argument("--panel-bed",   required=True, type=Path)
    ap.add_argument("--dictionary",  required=True, type=Path)
    ap.add_argument("--cytoband-bed", required=True, type=Path,
                    help="T2T-CHM13v2.0 cytoband BED (chrom start end band ...). "
                         "Off-panel breakpoint partners are characterized by band.")
    ap.add_argument("--sample",      required=True, type=str)
    ap.add_argument("--output",      required=True, type=Path)
    ap.add_argument("--version",     action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel_bed)
    dictionary = load_dictionary(args.dictionary)
    records = parse_vcf(args.vcf)
    cytobands = load_cytobands(args.cytoband_bed)
    rows = annotate(records, panel, dictionary, args.sample, cytobands)

    columns = [
        "sample", "sv_id", "sv_type", "filter",
        "chrom_a", "pos_a", "gene_a",
        "chrom_b", "pos_b", "gene_b",
        "gene_a_source", "gene_b_source",
        "known_mm_pair", "known_freq",
        "callers", "n_callers", "supp_vec", "support_reads",
    ]
    with open(args.output, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(c, "") for c in columns) + "\n")

    print(f"Annotated {len(rows)} on-panel SV records -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
