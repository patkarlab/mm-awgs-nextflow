#!/usr/bin/env python3
"""
build_ig_segments.py
====================

Emit a BED of immunoglobulin locus sub-regions: constant, J, D and V.

Why this exists
---------------
An IGH breakend's position within the locus carries mechanistic information
that its coordinate alone does not. Primary translocations in plasma cell
neoplasms arise from errors of class-switch recombination, which acts at the
switch regions immediately 5' of each constant-region gene. They therefore
break in the constant/switch region or near the J segments. The V array is
where V(D)J recombination and somatic hypermutation act physiologically, and
it is 121 near-identical segments over 940 kb, which is where reads mismap.

On this cohort, 82% of IGH breakends below the D-to-V transition carry a
dictionary name and 1.5% above it do. That is a mechanism, not a threshold,
and it should be read off the annotation rather than hardcoded as a
coordinate: the boundary moves between assemblies, and asserting it would
make the rule wrong on the next reference.

The segments are already in the RefSeq GFF as C_gene_segment,
J_gene_segment, D_gene_segment and V_gene_segment. This script collapses
each class to its extent within each Ig locus and writes it as a BED.

The constant-region extent is deliberately taken from the first base of the
locus rather than from the first C segment, because the switch regions
themselves are not annotated as features and sit 5' of the C segments they
serve. A breakend in Smu is in the switch region, not upstream of anything.

Nothing here is disease-specific and no coordinate is written by hand.

Usage
-----
  build_ig_segments.py \\
      --gff  GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz \\
      --panel-bed-chr assets/aWGS_PCN_v7_t2t_chr.bed \\
      --panel-bed-nc  assets/aWGS_PCN_v7_t2t_NC.bed \\
      --output assets/ig_segments_t2t.bed

Output columns: chrom, start, end, <LOCUS>_<CLASS>
e.g.  chr14  99839468  100127783  IGH_constant
"""

from __future__ import annotations

import argparse
import gzip
import re
import sys
from pathlib import Path

__version__ = "0.1.0"

NAME_RE = re.compile(r"(?:^|;)Name=([^;]+)")
# Segment features carry gene= and standard_name= rather than Name=, which
# only the parent gene feature has. Both are tried.
GENE_RE = re.compile(r"(?:^|;)gene=([^;]+)")

# Locus symbol -> the GFF gene feature that defines its extent.
LOCI = ("IGH", "IGK", "IGL")

CLASSES = {
    "C_gene_segment": "constant",
    "J_gene_segment": "J",
    "D_gene_segment": "D",
    "V_gene_segment": "V",
}


def open_maybe_gz(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def derive_rename(chr_bed: Path, nc_bed: Path) -> dict:
    """NC_ accession -> chr, by pairing the two panel BEDs line by line."""
    mapping = {}
    with open(chr_bed) as fc, open(nc_bed) as fn:
        for i, (lc, ln) in enumerate(zip(fc, fn), 1):
            c = lc.rstrip("\n").split("\t")
            n = ln.rstrip("\n").split("\t")
            if len(c) < 3 or len(n) < 3:
                continue
            if c[1:3] != n[1:3]:
                sys.exit(f"ERROR: panel BEDs disagree at line {i}")
            mapping[n[0]] = c[0]
    if not mapping:
        sys.exit("ERROR: no contig mapping derived")
    return mapping


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gff", required=True, type=Path)
    ap.add_argument("--panel-bed-chr", required=True, type=Path)
    ap.add_argument("--panel-bed-nc", required=True, type=Path)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    rename = derive_rename(args.panel_bed_chr, args.panel_bed_nc)

    # Pass 1: locus extents, from the gene feature carrying the locus symbol.
    extent = {}
    seg = {}          # (locus, class) -> [min, max]
    with open_maybe_gz(args.gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom = rename.get(f[0])
            if chrom is None:
                continue
            start, end = int(f[3]) - 1, int(f[4])
            if f[2] == "gene":
                m = NAME_RE.search(f[8])
                if m and m.group(1).upper() in LOCI:
                    extent[m.group(1).upper()] = (chrom, start, end)
                continue
            cls = CLASSES.get(f[2])
            if cls is None:
                continue
            m = NAME_RE.search(f[8]) or GENE_RE.search(f[8])
            sym = (m.group(1) if m else "").upper()
            locus = next((l for l in LOCI if sym.startswith(l)), None)
            if locus is None:
                continue
            key = (locus, cls)
            if key not in seg:
                seg[key] = [chrom, start, end]
            else:
                if seg[key][0] != chrom:
                    continue
                seg[key][1] = min(seg[key][1], start)
                seg[key][2] = max(seg[key][2], end)

    rows = []
    for locus in LOCI:
        if locus not in extent:
            sys.stderr.write(f"WARNING: no {locus} locus feature in the GFF\n")
            continue
        lchrom, lstart, lend = extent[locus]
        got = {}
        for cls in ("constant", "J", "D", "V"):
            v = seg.get((locus, cls))
            # Segments must lie on the locus's own chromosome and inside its
            # extent. Orphons carry the locus symbol but sit elsewhere: an
            # IGLV orphon on chr8 matches IGL by prefix and would otherwise
            # be emitted as the IGL V array.
            if not v or v[0] != lchrom or v[2] <= lstart or v[1] >= lend:
                continue
            got[cls] = (max(v[1], lstart), min(v[2], lend))

        # The constant class is extended from the locus start to the first
        # recombination segment. Switch regions are not annotated as
        # features and lie between the C genes and JH, so the gap would
        # otherwise be unclassified -- and on this cohort that gap is where
        # the t(11;14) and t(4;14) breakpoints sit. Class-switch
        # recombination acts there, which is why primary translocations do.
        vdj_start = min((s for c, (s, e) in got.items() if c != "constant"),
                        default=None)
        if "constant" in got or vdj_start is not None:
            c_end = vdj_start if vdj_start is not None else got["constant"][1]
            rows.append((lchrom, lstart, c_end, f"{locus}_C_switch"))
        for cls in ("J", "D", "V"):
            if cls in got:
                s, e = got[cls]
                rows.append((lchrom, s, e, f"{locus}_{cls}"))

    rows.sort(key=lambda r: (r[0], r[1]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        for chrom, s, e, name in rows:
            out.write(f"{chrom}\t{s}\t{e}\t{name}\n")

    for chrom, s, e, name in rows:
        sys.stderr.write(f"  {name:14} {chrom}:{s:,}-{e:,}  ({e - s:,} bp)\n")
    sys.stderr.write(f"wrote {len(rows)} segment class(es) -> {args.output}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
