#!/usr/bin/env python3
"""
build_mm_gene_model.py
======================

Emit one BED feature per gene named by the panel: bare gene bodies, no flank,
no merging.

Why this exists
---------------
The panel BED serves capture. Its intervals are flanked to cover breakpoint
cluster regions and overlapping ones are merged, which is correct for
enrichment and wrong for naming a breakpoint. Two consequences:

  Compound labels. A t(4;14) breakend lands in a window labelled
  "FGFR3/NSD2", so the annotated table cannot say which of the two it is,
  and IGH::FGFR3 cannot be distinguished from IGH::NSD2. The same applies to
  WWOX/MAF, TP53+TNFSF12, FCRL5/FCRL4 and IGL/IGLL5. In this reference FGFR3
  and NSD2 are 62 kb apart and WWOX and MAF do not overlap at all, so the
  distinction is available and is simply being discarded.

  No distance. 62% of the v7 panel is flank rather than gene body, by design:
  the BCR windows are what make Ig translocations detectable. But it means a
  breakend labelled MYC may be anywhere in a 5 Mb window, and the table gives
  no way to tell a breakend in the gene from one 2 Mb away.

This model is for annotation only. It never gates panel membership, which
stays with the panel BED. It answers one question: which gene is this
coordinate in, and how far is it from the gene it was named after.

Where the coordinates come from
-------------------------------
The RefSeq GFF for the same assembly, not from this script. Nothing here
carries a coordinate. The symbols come from column 4 of the panel BED, split
on "/" and "+", so the model tracks the panel automatically and a panel
revision needs no edit here.

Ig loci are included. IGH, IGK and IGL are annotated in RefSeq as locus
features of biotype "other" with real spans, so features are taken by symbol
match regardless of biotype. IGLL5 sits inside the IGL locus; the annotator
resolves a coordinate to the tightest containing feature, so a breakend in
IGLL5 is named IGLL5 and one elsewhere in the locus is named IGL.

Contig naming
-------------
The NC_ to chr mapping is derived by pairing the two panel BEDs line by line,
which the panel README states carry identical coordinates and differ only in
the contig column. Deriving it means there is no third file to keep in step,
and a mismatch between the two BEDs is detected here rather than silently
producing a model on the wrong contigs.

Usage
-----
  build_mm_gene_model.py \\
      --panel-bed-chr assets/aWGS_PCN_v7_t2t_chr.bed \\
      --panel-bed-nc  assets/aWGS_PCN_v7_t2t_NC.bed \\
      --gff           GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz \\
      --output        assets/aWGS_PCN_v7_gene_model_t2t.bed

Add symbols the panel does not name with --extra-genes SYM[,SYM...], for a
partner that must be nameable but is not a capture target.
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

__version__ = "0.1.0"

NAME_RE = re.compile(r"(?:^|;)Name=([^;]+)")
BIOTYPE_RE = re.compile(r"(?:^|;)gene_biotype=([^;]+)")


def open_maybe_gz(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def derive_rename(chr_bed: Path, nc_bed: Path) -> dict:
    """NC_ accession -> chr name, from the two paired panel BEDs.

    The panel README states the files carry identical coordinates and differ
    only in the contig column, so pairing them line by line is enough. Any
    disagreement in the coordinate columns is fatal: it would mean the two
    BEDs have drifted, and a gene model built on the wrong contig naming
    would produce plausible coordinates on the wrong chromosomes.
    """
    mapping = {}
    with open(chr_bed) as fc, open(nc_bed) as fn:
        for i, (lc, ln) in enumerate(zip(fc, fn), 1):
            c = lc.rstrip("\n").split("\t")
            n = ln.rstrip("\n").split("\t")
            if len(c) < 3 or len(n) < 3:
                continue
            if c[1:3] != n[1:3]:
                sys.exit(f"ERROR: {chr_bed.name} and {nc_bed.name} disagree at "
                         f"line {i}: {c[1:3]} vs {n[1:3]}. The two BEDs must "
                         f"carry identical coordinates.")
            prev = mapping.get(n[0])
            if prev is not None and prev != c[0]:
                sys.exit(f"ERROR: {n[0]} maps to both {prev} and {c[0]}")
            mapping[n[0]] = c[0]
    if not mapping:
        sys.exit("ERROR: no contig mapping derived from the panel BEDs")
    return mapping


def panel_symbols(chr_bed: Path) -> set:
    """Gene symbols named by the panel, from column 4.

    Compound labels are split on "/" and "+", the two joiners the panel uses:
    FGFR3/NSD2, WWOX/MAF, FCRL5/FCRL4, IGL/IGLL5, TP53+TNFSF12.
    """
    out = set()
    with open(chr_bed) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            for tok in parts[3].replace("+", "/").split("/"):
                tok = tok.strip()
                if tok:
                    out.add(tok.upper())
    if not out:
        sys.exit(f"ERROR: no symbols in column 4 of {chr_bed}")
    return out


def collect_features(gff: Path, wanted: set, rename: dict):
    """Every GFF gene feature whose Name is a wanted symbol.

    Biotype is recorded but never filtered on. IGH, IGK and IGL are biotype
    "other" and are exactly the features this panel exists to anchor against.
    """
    found = {}
    with open_maybe_gz(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom = rename.get(f[0])
            if chrom is None:
                continue
            m = NAME_RE.search(f[8])
            if not m:
                continue
            sym = m.group(1)
            if sym.upper() not in wanted:
                continue
            b = BIOTYPE_RE.search(f[8])
            found.setdefault(sym.upper(), []).append(
                (chrom, int(f[3]) - 1, int(f[4]), sym, b.group(1) if b else ""))
    return found


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--panel-bed-chr", required=True, type=Path)
    ap.add_argument("--panel-bed-nc", required=True, type=Path)
    ap.add_argument("--gff", required=True, type=Path,
                    help="RefSeq GFF for the same assembly; .gz accepted.")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--extra-genes", default="",
                    help="Comma-separated symbols to include beyond those the "
                         "panel names.")
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    rename = derive_rename(args.panel_bed_chr, args.panel_bed_nc)
    sys.stderr.write(f"contig mapping: {len(rename)} accessions "
                     f"derived from the paired panel BEDs\n")

    wanted = panel_symbols(args.panel_bed_chr)
    extra = {s.strip().upper() for s in args.extra_genes.split(",") if s.strip()}
    wanted |= extra
    sys.stderr.write(f"symbols wanted: {len(wanted)} "
                     f"({len(extra)} beyond the panel)\n")

    found = collect_features(args.gff, wanted, rename)

    rows = []
    duplicated = []
    for sym in sorted(found):
        feats = found[sym]
        if len(feats) > 1:
            # More than one record for a symbol. Keep them all rather than
            # choosing: the annotator resolves a coordinate to the tightest
            # containing feature and measures distance to the nearest record
            # of that name, so several records are handled correctly and
            # picking one here would discard a real locus.
            duplicated.append(f"{sym} ({len(feats)})")
        rows.extend(feats)

    missing = sorted(wanted - set(found))
    if missing:
        sys.stderr.write(
            "WARNING: no GFF gene feature for: " + ", ".join(missing) + "\n"
            "  These stay nameable through the panel interval label, but a\n"
            "  breakend in one cannot be distinguished from its neighbours in\n"
            "  the same window and carries no distance.\n")
    if duplicated:
        sys.stderr.write("symbols with multiple records, all kept: "
                         + ", ".join(duplicated) + "\n")

    def sort_key(r):
        c = r[0][3:]
        return (int(c) if c.isdigit() else {"X": 98, "Y": 99}.get(c, 100), r[1])

    rows.sort(key=sort_key)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        for chrom, start, end, sym, _bio in rows:
            out.write(f"{chrom}\t{start}\t{end}\t{sym}\n")

    span = sum(e - s for _c, s, e, _n, _b in rows)
    sys.stderr.write(
        f"wrote {len(rows)} features for {len(found)} symbols, "
        f"{span:,} bp -> {args.output}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
