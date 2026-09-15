#!/usr/bin/env python3
"""
ichorcna_to_mmplot.py

Convert ichorCNA output into the two positional TSVs mmplot.py expects, so the
chromosome figures keep the lab's house style while the copy number behind
them comes from ichorCNA rather than mmkaryo.

What mmplot needs
-----------------
mmplot takes two positional arguments before the output path:

  bins      columns chromosome, start, log2_ratio
  segments  columns chromosome, start, end, log2_ratio, and optionally theta

Both are read by name from the header, so extra columns are ignored and a
missing optional one is treated as absent rather than an error.

Mapping
-------
  .cna.seg  -> bins.tsv       one row per 1 Mb bin
  .seg.txt  -> segments.tsv   one row per segment

ichorCNA writes the sample id as a prefix on the .cna.seg column names
(11F20265231.logR), so columns are matched by suffix. Bins whose logR is NA
are written through as NA: mmplot skips them, and dropping them here instead
would silently close the gaps where the panel and centromere masks are, which
is exactly where a reader should see nothing rather than an interpolated line.

Allele fraction
---------------
No theta column is written. ichorCNA is depth-only and cannot separate 2:0
from 1:1, and under the current scope allelic state is measured at eight
panel windows rather than genome-wide, so there is no per-bin allele fraction
to supply. mmplot tolerates the absence: call_class falls back to
DEL / balanced het / GAIN from total copy number alone, and minor allele copy
number comes out NaN so its trace draws nothing.

That leaves an empty BAF panel in the figure, which reads as measured and
found nothing rather than not measured. Removing the panel is a change to
mmplot itself and is deliberately not done here.

Dependencies: standard library only.
"""

import argparse
import csv
import os
import sys

MAIN_CHROMS = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]


def column_by_suffix(fieldnames, *suffixes):
    """First column matching a suffix, allowing ichorCNA's sample-id prefix."""
    for suffix in suffixes:
        for name in fieldnames or []:
            if name == suffix or name.endswith("." + suffix):
                return name
    return None


def convert_bins(cna_path, out_path):
    """.cna.seg to the bins table, preserving NA rows."""
    written = na = 0
    with open(cna_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        logr_col = column_by_suffix(reader.fieldnames, "logR")
        cn_col = column_by_suffix(reader.fieldnames,
                                  "Corrected_Copy_Number", "copy.number")
        if logr_col is None:
            sys.exit("ERROR: no logR column in %s" % cna_path)

        with open(out_path, "w", encoding="utf-8") as out:
            out.write("chromosome\tstart\tend\tlog2_ratio\tcopy_number\n")
            for row in reader:
                chrom = row.get("chr") or row.get("chrom") or ""
                if chrom not in MAIN_CHROMS:
                    continue
                try:
                    start, end = int(row["start"]), int(row["end"])
                except (KeyError, TypeError, ValueError):
                    continue
                raw = (row.get(logr_col) or "").strip()
                try:
                    value = float(raw)
                    if value != value:
                        raise ValueError
                    log2 = "%.4f" % value
                except ValueError:
                    log2, na = "NA", na + 1
                copies = (row.get(cn_col) or "").strip() if cn_col else ""
                try:
                    copies = str(int(float(copies)))
                except ValueError:
                    copies = "NA"
                out.write("%s\t%d\t%d\t%s\t%s\n" % (chrom, start, end, log2, copies))
                written += 1
    return written, na


def convert_segments(seg_path, out_path):
    """.seg.txt to the segments table."""
    written = 0
    with open(seg_path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        logr_col = column_by_suffix(reader.fieldnames,
                                    "seg.median.logR", "logR")
        cn_col = column_by_suffix(reader.fieldnames,
                                  "Corrected_Copy_Number", "copy.number")
        call_col = column_by_suffix(reader.fieldnames,
                                    "Corrected_Call", "call")
        bins_col = column_by_suffix(reader.fieldnames, "num.mark")
        if logr_col is None:
            sys.exit("ERROR: no segment log2 column in %s" % seg_path)

        with open(out_path, "w", encoding="utf-8") as out:
            out.write("chromosome\tstart\tend\tlog2_ratio\tcopy_number\t"
                      "call\tn_bins\n")
            for row in reader:
                chrom = row.get("chrom") or row.get("chr") or ""
                if chrom not in MAIN_CHROMS:
                    continue
                try:
                    start, end = int(row["start"]), int(row["end"])
                    log2 = float(row[logr_col])
                except (KeyError, TypeError, ValueError):
                    continue
                if log2 != log2:
                    continue
                copies = (row.get(cn_col) or "").strip() if cn_col else ""
                try:
                    copies = str(int(float(copies)))
                except ValueError:
                    copies = "NA"
                call = ((row.get(call_col) or "").strip() or "NA") if call_col else "NA"
                n_bins = ((row.get(bins_col) or "").strip() or "0") if bins_col else "0"
                out.write("%s\t%d\t%d\t%.4f\t%s\t%s\t%s\n"
                          % (chrom, start, end, log2, copies, call, n_bins))
                written += 1
    return written


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cna", required=True,
                        help="ichorCNA <sample>.cna.seg")
    parser.add_argument("--seg", required=True,
                        help="ichorCNA <sample>.seg.txt")
    parser.add_argument("--out-prefix", required=True,
                        help="writes <prefix>.bins.tsv and <prefix>.segments.tsv")
    args = parser.parse_args()

    for path in (args.cna, args.seg):
        if not os.path.isfile(path):
            sys.exit("ERROR: %s not found" % path)

    bins_path = args.out_prefix + ".bins.tsv"
    segments_path = args.out_prefix + ".segments.tsv"

    n_bins, n_na = convert_bins(args.cna, bins_path)
    n_segments = convert_segments(args.seg, segments_path)

    if not n_bins:
        sys.exit("ERROR: no usable bins written from %s" % args.cna)
    if not n_segments:
        sys.exit("ERROR: no usable segments written from %s" % args.seg)

    print("bins      : %d (%d NA, kept as NA so the masked regions stay blank)"
          % (n_bins, n_na))
    print("segments  : %d" % n_segments)
    print("written   : %s" % bins_path)
    print("            %s" % segments_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
