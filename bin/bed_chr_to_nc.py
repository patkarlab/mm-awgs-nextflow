#!/usr/bin/env python3
"""
bed_chr_to_nc.py

Write an NC_-named copy of a chr-named BED, taking the chr -> NC_ accession
mapping from the RefSeq GFF of the same assembly rather than from a
hand-typed table.

Why: build_mm_gene_model.py derives contig naming by pairing a chr-named and
an NC_-named panel BED line by line. The T2T panel ships both; the hg38 panel
ships only the chr-named file. Typing GRCh38.p14 accession versions by hand
(NC_000001.11, NC_000014.9, ...) is where a silent wrong version would come
from, so the mapping is read from the GFF `region` records, which carry
`chromosome=N` on the primary-assembly contigs.

Only primary-assembly chromosomes are mapped (NC_ accessions). Alt, patch and
unlocalized contigs (NT_/NW_) never appear in a panel BED and are ignored.

Usage:
  bed_chr_to_nc.py --gff <RefSeq .gff or .gff.gz> --bed <chr BED> --output <NC BED>
"""

import argparse
import gzip
import sys


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def chr_to_nc_from_gff(gff_path):
    """Return {'chr1': 'NC_000001.11', ..., 'chrX': ..., 'chrY': ..., 'chrM': ...}."""
    mapping = {}
    with open_text(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "region":
                continue
            seqid = fields[0]
            if not seqid.startswith("NC_"):
                continue
            attrs = dict(kv.split("=", 1) for kv in fields[8].split(";") if "=" in kv)
            chrom = attrs.get("chromosome")
            if chrom is None:
                continue
            # Mitochondrion is annotated chromosome=MT; UCSC naming is chrM.
            name = "chrM" if chrom == "MT" else "chr" + chrom
            mapping.setdefault(name, seqid)
            # All 25 primary contigs are collected long before the end of a
            # 78 MB file; stop once they are in hand.
            if len(mapping) >= 25:
                break
    return mapping


def main():
    ap = argparse.ArgumentParser(
        description="NC_-named copy of a chr-named BED via RefSeq GFF region records.")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--bed", required=True, help="chr-named input BED")
    ap.add_argument("--output", required=True, help="NC_-named output BED")
    args = ap.parse_args()

    mapping = chr_to_nc_from_gff(args.gff)
    if len(mapping) < 24:
        sys.exit("ERROR: only %d primary contigs found in %s; expected 24-25."
                 % (len(mapping), args.gff))

    n = 0
    unmapped = set()
    with open(args.bed) as fin, open(args.output, "w") as fout:
        for line in fin:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            nc = mapping.get(fields[0])
            if nc is None:
                unmapped.add(fields[0])
                continue
            fields[0] = nc
            fout.write("\t".join(fields) + "\n")
            n += 1
    if unmapped:
        sys.exit("ERROR: contigs with no NC_ accession in the GFF: %s"
                 % ", ".join(sorted(unmapped)))

    print("Mapping (from GFF region records):")
    for k in sorted(mapping, key=lambda c: (len(c), c)):
        print("  %-5s %s" % (k, mapping[k]))
    print("Wrote %d regions to %s" % (n, args.output))


if __name__ == "__main__":
    main()
