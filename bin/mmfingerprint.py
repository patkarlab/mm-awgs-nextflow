#!/usr/bin/env python3
"""
mmfingerprint.py

Is this BAM one person? Allele fractions at common SNPs on the panel, read
from on-target reads, answer that independently of the sex field on the
sample sheet.

A single individual, at a site of adequate depth, shows an alternate allele
fraction near 0, near 0.5 or near 1. Copy-number change moves heterozygous
sites off 0.5 in one direction per region (a trisomy at clonality c puts
them at (1+c)/(2+c) and its complement), but a diploid region stays at 0.5.
A two-person mixture at fraction f does something no single genome does: at
every site where the two genotypes differ it puts the allele fraction at f/2,
f, 1-f or 1-f/2, uniformly across the genome, diploid regions included, and
the count of such sites is large because two unrelated people differ at a
third of common SNPs.

So the statistic reported per window and genome-wide is the fraction of
usable sites with an intermediate allele fraction, in the bands
--low..--mid-low and --mid-high..--high (defaults 0.08-0.35 and 0.65-0.92),
together with the observed heterozygous fraction against the population
expectation. A clean sample has a few percent intermediate sites (sequencing
error and subclonal CNA); a mixture has tens of percent. Sex is read from
chrX heterozygosity and chrY depth relative to autosomes. Loss of Y in a male
tumour drops chrY depth and leaves everything else alone; a mixture does not.

Interpretation is left to the reader: the script prints the numbers and a
one-line summary, and compares against a second BAM if given, so the
suspect sample can be read beside a clean one from the same run.

Dependencies: pysam (conda activate awgs_sv); imports genebaf.py from the
same directory for its loaders.

Usage:
  mmfingerprint.py <bam> <panel_bed> <sites_vcf> --sample ID [--compare <bam2>]
      [--af-field SAS_AF] [--min-depth 15] [--min-alt 2]
      [--out <prefix>]     writes <prefix>.fingerprint.tsv and .json
"""

import argparse
import json
import os
import statistics
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
try:
    import genebaf
except ImportError:
    sys.exit("ERROR: genebaf.py must sit beside this script")
try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is required (conda activate awgs_sv)")


def per_site_fractions(bam_path, windows, sites, min_depth, min_alt, min_mapq, min_baseq,
                       max_other_fraction):
    """[(chrom, pos, af, depth, is_het, expected_het)] over all windows."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    out = []
    for window in windows:
        chrom, start, end, name = window
        site_list = sites.get(window, [])
        if not site_list:
            continue
        counts = genebaf.tally_window(bam, chrom, start, end, site_list, min_mapq, min_baseq)
        for pos, ref, alt, af_pop in site_list:
            r, a, o = counts[pos]
            n = r + a
            if n < min_depth or o > max_other_fraction * (n + o):
                continue
            af = a / n
            is_het = (r >= min_alt and a >= min_alt)
            out.append((name, chrom, pos, af, n, is_het, 2 * af_pop * (1 - af_pop)))
    bam.close()
    return out


def chrom_depth(bam_path, chroms, n_probe=200, probe_len=100000):
    """Mean depth over evenly spaced probes on each chromosome, so chrY and
    chrX can be read against autosomes without a panel window on them."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    lengths = dict(zip(bam.references, bam.lengths))
    depth = {}
    for chrom in chroms:
        if chrom not in lengths:
            depth[chrom] = None
            continue
        L = lengths[chrom]
        step = max(probe_len, L // n_probe)
        total_bases = 0
        total_len = 0
        for start in range(int(L * 0.05), int(L * 0.95), step):
            end = min(start + probe_len, L)
            cov = bam.count_coverage(chrom, start, end, quality_threshold=0, read_callback="nofilter")
            total_bases += sum(sum(c) for c in cov)
            total_len += end - start
        depth[chrom] = (total_bases / total_len) if total_len else None
    bam.close()
    return depth


def summarise(rows, bands):
    low, mid_low, mid_high, high = bands
    usable = len(rows)
    inter = sum(1 for r in rows if (low <= r[3] <= mid_low) or (mid_high <= r[3] <= high))
    het = sum(1 for r in rows if r[5])
    exp = statistics.mean(r[6] for r in rows) if rows else 0.0
    obs = het / usable if usable else 0.0
    return {
        "usable_sites": usable,
        "intermediate_sites": inter,
        "intermediate_fraction": round(inter / usable, 4) if usable else None,
        "observed_het_fraction": round(obs, 4) if usable else None,
        "expected_het_fraction": round(exp, 4) if usable else None,
        "het_ratio": round(obs / exp, 3) if (usable and exp) else None,
    }


def histogram(rows, nbins=20):
    counts = [0] * nbins
    for r in rows:
        counts[min(nbins - 1, int(r[3] * nbins))] += 1
    return counts


def fingerprint(bam, windows, sites, args):
    rows = per_site_fractions(bam, windows, sites, args.min_depth, args.min_alt,
                              args.min_mapq, args.min_baseq, args.max_other_fraction)
    bands = (args.low, args.mid_low, args.mid_high, args.high)
    autosomal = [r for r in rows if r[1] not in ("chrX", "chrY")]
    x_rows = [r for r in rows if r[1] == "chrX"]
    per_window = {}
    for name in sorted({r[0] for r in rows}):
        per_window[name] = summarise([r for r in rows if r[0] == name], bands)
    depth = chrom_depth(bam, ["chr1", "chr2", "chr3", "chrX", "chrY"])
    auto_depth = [d for c, d in depth.items() if c in ("chr1", "chr2", "chr3") and d]
    auto_mean = statistics.mean(auto_depth) if auto_depth else None
    x_ratio = (depth["chrX"] / auto_mean) if (auto_mean and depth.get("chrX")) else None
    y_ratio = (depth["chrY"] / auto_mean) if (auto_mean and depth.get("chrY") is not None) else None
    x_summary = summarise(x_rows, bands) if x_rows else None

    # Sex from depth: chrX near 1.0 of autosomal and chrY near 0 is XX;
    # chrX near 0.5 and chrY near 0.5 is XY; chrX near 0.5 and chrY near 0
    # is XY with loss of Y (or a male whose Y is not covered).
    if x_ratio is None:
        sex_call = "unknown"
    elif x_ratio > 0.75:
        sex_call = "XX"
    elif y_ratio is not None and y_ratio > 0.2:
        sex_call = "XY"
    else:
        sex_call = "XY_loss_of_Y_or_unread_Y"

    return {
        "autosomal": summarise(autosomal, bands),
        "autosomal_af_histogram_20bins": histogram(autosomal),
        "chrX": x_summary,
        "depth_ratio_chrX_over_autosome": round(x_ratio, 3) if x_ratio is not None else None,
        "depth_ratio_chrY_over_autosome": round(y_ratio, 3) if y_ratio is not None else None,
        "sex_from_depth": sex_call,
        "per_window": per_window,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bam")
    ap.add_argument("panel_bed")
    ap.add_argument("sites_vcf")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--compare", default=None, help="a second BAM to report beside the first")
    ap.add_argument("--compare-sample", default="comparison")
    ap.add_argument("--out", default=None, help="output prefix [sample]")
    ap.add_argument("--af-field", default="SAS_AF")
    ap.add_argument("--min-af", type=float, default=0.25)
    ap.add_argument("--max-af", type=float, default=0.75)
    ap.add_argument("--min-depth", type=int, default=15)
    ap.add_argument("--min-alt", type=int, default=2)
    ap.add_argument("--min-mapq", type=int, default=10)
    ap.add_argument("--min-baseq", type=int, default=10)
    ap.add_argument("--max-other-fraction", type=float, default=0.2)
    ap.add_argument("--low", type=float, default=0.08)
    ap.add_argument("--mid-low", type=float, default=0.35)
    ap.add_argument("--mid-high", type=float, default=0.65)
    ap.add_argument("--high", type=float, default=0.92)
    args = ap.parse_args()

    windows = genebaf.load_bed(args.panel_bed)
    sites = genebaf.load_sites_in_windows(args.sites_vcf, windows, args.af_field,
                                          args.min_af, args.max_af)
    sys.stderr.write(f"{len(windows)} windows, {sum(len(v) for v in sites.values())} common SNPs\n")

    result = {"sample": args.sample, "bam": args.bam,
              "parameters": {k: getattr(args, k) for k in
                             ("min_depth", "min_alt", "low", "mid_low", "mid_high", "high",
                              "af_field", "min_af", "max_af")}}
    result["primary"] = fingerprint(args.bam, windows, sites, args)
    if args.compare:
        result["comparison"] = fingerprint(args.compare, windows, sites, args)
        result["comparison_sample"] = args.compare_sample

    prefix = args.out or args.sample
    with open(prefix + ".fingerprint.json", "w", encoding="utf-8") as fh:
        json.dump(result, fh, indent=2)
    with open(prefix + ".fingerprint.tsv", "w", encoding="utf-8") as fh:
        fh.write("sample\twindow\tusable_sites\tintermediate_sites\tintermediate_fraction\t"
                 "observed_het_fraction\texpected_het_fraction\thet_ratio\n")
        for label, block in ((args.sample, result["primary"]),
                             (args.compare_sample, result.get("comparison"))):
            if not block:
                continue
            for name, s in [("AUTOSOMAL", block["autosomal"])] + sorted(block["per_window"].items()):
                fh.write("\t".join(str(x) for x in (
                    label, name, s["usable_sites"], s["intermediate_sites"],
                    s["intermediate_fraction"], s["observed_het_fraction"],
                    s["expected_het_fraction"], s["het_ratio"])) + "\n")

    def line(label, block):
        a = block["autosomal"]
        print(f"{label:16s} usable {a['usable_sites']:6d}  intermediate {a['intermediate_fraction']!s:>7}  "
              f"het obs/exp {a['observed_het_fraction']}/{a['expected_het_fraction']} "
              f"(ratio {a['het_ratio']})  chrX/auto {block['depth_ratio_chrX_over_autosome']}  "
              f"chrY/auto {block['depth_ratio_chrY_over_autosome']}  sex_from_depth {block['sex_from_depth']}")
        print(f"{'':16s} AF histogram (20 bins 0..1): {block['autosomal_af_histogram_20bins']}")

    line(args.sample, result["primary"])
    if args.compare:
        line(args.compare_sample, result["comparison"])
    print(f"written: {prefix}.fingerprint.tsv, {prefix}.fingerprint.json")


if __name__ == "__main__":
    main()
