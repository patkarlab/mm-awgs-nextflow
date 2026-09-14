#!/usr/bin/env python3
"""
genebaf.py

Allelic state per panel gene, from on-target reads.

Why on-target only
------------------
Off-target coverage in this cohort is 0.6 to 3x. At that depth most common SNPs
are seen in one or two reads, and a heterozygous site shows one allele by chance
about half the time. A likelihood that marginalises over genotypes handles the
algebra correctly but has almost no power to distinguish balanced from
imbalanced, and any excess of homozygous-looking sites drags the estimate toward
LOH. Measured on one sample, 1 Mb windows on a chromosome the depth caller
reports as diploid returned theta anywhere from 0.01 to 0.50 with no centre; 5 Mb
windows at 300 to 800 sites did no better. Off-target BAF can see a complete LOH
and nothing else.

On-target coverage is 19 to 21x on properly sequenced samples. At that depth a
het is unambiguous, and the unfolded allele fraction at hets reads the allelic
ratio directly. On one sample the deep hets across three chr1 panel windows were
bimodal at 0.2-0.3 and 0.6-0.7 with a dip at 0.4 — a 2:1 ratio read straight
from the reads, matching the depth caller's 1p loss at 77 percent cell fraction.

So: depth genome-wide from ichorCNA on off-target reads; allelic state per gene
from on-target reads here. Each method is used where the data support it.

What is measured
----------------
For each panel window the script pileups every common SNP from the supplied
list, keeps sites with at least --min-depth reads, and classifies each as
homozygous (fewer than --min-alt reads on either allele) or heterozygous. At
hets it records the unfolded alt fraction. From the het set:

  n_het            hets with enough depth to use
  median_dev       median |af - 0.5|, the folded distance from balance
  bimodal_gap      fraction of hets in the 0.40-0.60 band; low means the
                   distribution has split into two modes
  allelic_state    balanced / imbalanced / loh / insufficient

  LOH is called when almost no hets survive relative to the number of
  homozygous sites at usable depth: a deleted or copy-neutral-LOH region has
  the population's het sites but the sample shows one allele at all of them.
  A comparison against the expected het fraction from the population AF
  makes this a ratio rather than a raw count.

Depth is reported alongside so a reader can pair the allelic call with the
copy number from ichorkaryo: depth-neutral plus LOH is copy-neutral LOH;
depth-loss plus LOH is a deletion; depth-gain plus 2:1 is a trisomy.

Dependencies: pysam.
"""

import argparse
import gzip
import json
import statistics
import sys
from collections import defaultdict

try:
    import pysam
except ImportError:
    sys.exit("ERROR: pysam is required (conda activate awgs_sv)")

MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def canonical(chrom):
    """Return chr-prefixed name for a primary chromosome, else None."""
    name = chrom if chrom.startswith("chr") else "chr" + chrom
    if name == "chrMT":
        name = "chrM"
    return name if name in MAIN_CHROMS else None


def load_bed(path):
    """Panel windows as [(chrom, start, end, name), ...]."""
    windows = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            chrom = canonical(parts[0])
            if chrom is None or len(parts) < 3:
                continue
            name = parts[3] if len(parts) > 3 else f"{chrom}:{parts[1]}-{parts[2]}"
            windows.append((chrom, int(parts[1]), int(parts[2]), name))
    return windows


def load_sites_in_windows(vcf_path, windows, af_field, min_af, max_af):
    """{(chrom,start,end,name): [(pos, ref, alt, af), ...]} for sites inside windows."""
    by_chrom = defaultdict(list)
    for w in windows:
        by_chrom[w[0]].append(w)

    sites = defaultdict(list)
    opener = gzip.open if vcf_path.endswith(".gz") else open
    with opener(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 8)
            chrom = canonical(parts[0])
            if chrom is None or chrom not in by_chrom:
                continue
            pos = int(parts[1])
            ref, alt = parts[3], parts[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            af = None
            for kv in parts[7].split(";"):
                if kv.startswith(af_field + "="):
                    try:
                        af = float(kv.split("=", 1)[1].split(",")[0])
                    except ValueError:
                        af = None
                    break
            if af is None or not (min_af <= af <= max_af):
                continue
            for w in by_chrom[chrom]:
                if w[1] < pos <= w[2]:
                    sites[w].append((pos, ref, alt, af))
    return sites


def tally_window(bam, chrom, start, end, site_list, min_mapq, min_baseq):
    """Allele counts at every site in one window from a single pass.

    Returns {pos: [ref_count, alt_count, other_count]} keyed by 1-based
    position. Each read is visited once; its aligned pairs are walked and any
    reference position that is a site of interest is counted. This replaces
    one pileup() iterator per site with one fetch() per window, which is the
    difference between minutes and seconds per sample.
    """
    wanted = {pos: (ref, alt) for pos, ref, alt, _af in site_list}
    counts = {pos: [0, 0, 0] for pos in wanted}
    if not wanted:
        return counts
    lo, hi = min(wanted), max(wanted)

    for read in bam.fetch(chrom, lo - 1, hi):
        if read.is_unmapped or read.is_secondary or read.is_supplementary \
                or read.is_duplicate or read.is_qcfail:
            continue
        if read.mapping_quality < min_mapq:
            continue
        seq = read.query_sequence
        if seq is None:
            continue
        quals = read.query_qualities
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            pos1 = rpos + 1
            if pos1 not in wanted:
                continue
            if quals is not None and quals[qpos] < min_baseq:
                continue
            base = seq[qpos].upper()
            ref, alt = wanted[pos1]
            entry = counts[pos1]
            if base == ref:
                entry[0] += 1
            elif base == alt:
                entry[1] += 1
            else:
                entry[2] += 1
    return counts


def window_depth(bam, chrom, start, end):
    """Mean depth across a window from samtools-style coverage."""
    try:
        cov = bam.count_coverage(chrom, start, end, quality_threshold=0)
    except ValueError:
        return float("nan")
    total = sum(sum(c) for c in cov)
    length = end - start
    return total / length if length else float("nan")


def classify(n_het, n_hom_usable, expected_het_fraction, devs, in_band,
             min_hets, loh_ratio):
    """Allelic state from the het allele-fraction distribution."""
    usable = n_het + n_hom_usable
    if usable == 0:
        return "insufficient", None
    observed_het_fraction = n_het / usable
    # A region where hets have collapsed to one allele shows far fewer hets
    # than the population predicts. loh_ratio is how far below expectation
    # counts as loss of heterozygosity.
    if expected_het_fraction > 0 and \
            observed_het_fraction < expected_het_fraction * loh_ratio and usable >= min_hets:
        return "loh", round(observed_het_fraction / expected_het_fraction, 3)
    if n_het < min_hets:
        return "insufficient", None
    median_dev = statistics.median(devs)
    # Balanced hets scatter around 0.5 with a spread set by depth; at 20x the
    # binomial SD is about 0.11, so a median deviation past 0.15 combined with
    # most hets outside the central band is a split distribution, not noise.
    if median_dev >= 0.15 and in_band < 0.35:
        return "imbalanced", round(median_dev, 3)
    return "balanced", round(median_dev, 3)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", help="hg38 BAM, indexed")
    parser.add_argument("panel_bed", help="panel BED with gene names in column 4")
    parser.add_argument("sites_vcf", help="common biallelic SNP VCF with a population AF")
    parser.add_argument("out_prefix")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--af-field", default="SAS_AF")
    parser.add_argument("--min-af", type=float, default=0.25)
    parser.add_argument("--max-af", type=float, default=0.75)
    parser.add_argument("--min-depth", type=int, default=10,
                        help="reads required at a site before it is used [10]")
    parser.add_argument("--min-alt", type=int, default=3,
                        help="reads on the minor allele required to call a site het [3]")
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-baseq", type=int, default=10)
    parser.add_argument("--min-hets", type=int, default=8,
                        help="hets required for a balanced/imbalanced call [8]")
    parser.add_argument("--loh-ratio", type=float, default=0.25,
                        help="observed/expected het fraction below which LOH is "
                             "called [0.25]")
    parser.add_argument("--max-other-fraction", type=float, default=0.2,
                        help="drop sites where third-allele reads exceed this "
                             "fraction of depth [0.2]")
    args = parser.parse_args()

    windows = load_bed(args.panel_bed)
    if not windows:
        sys.exit(f"ERROR: no usable windows in {args.panel_bed}")
    sys.stderr.write(f"{len(windows)} panel windows\n")

    sites = load_sites_in_windows(args.sites_vcf, windows, args.af_field,
                                  args.min_af, args.max_af)
    total_sites = sum(len(v) for v in sites.values())
    sys.stderr.write(f"{total_sites} common SNPs inside panel windows\n")

    bam = pysam.AlignmentFile(args.bam, "rb")
    rows = []
    per_site_rows = []

    for window in windows:
        chrom, start, end, name = window
        depth = window_depth(bam, chrom, start, end)
        devs, afs = [], []
        n_het = n_hom = n_lowdepth = n_other = 0
        expected_het = []

        site_list = sites.get(window, [])
        counts = tally_window(bam, chrom, start, end, site_list,
                              args.min_mapq, args.min_baseq)
        for pos, ref, alt, af in site_list:
            r, a, o = counts[pos]
            n = r + a
            if n < args.min_depth:
                n_lowdepth += 1
                continue
            if o > args.max_other_fraction * (n + o):
                n_other += 1
                continue
            expected_het.append(2 * af * (1 - af))
            if r >= args.min_alt and a >= args.min_alt:
                n_het += 1
                frac = a / n
                afs.append(frac)
                devs.append(abs(frac - 0.5))
                per_site_rows.append((name, chrom, pos, r, a, round(frac, 3)))
            else:
                n_hom += 1

        expected_het_fraction = (statistics.mean(expected_het) if expected_het else 0.0)
        in_band = (sum(1 for f in afs if 0.40 <= f <= 0.60) / n_het) if n_het else 0.0
        state, stat = classify(n_het, n_hom, expected_het_fraction, devs, in_band,
                               args.min_hets, args.loh_ratio)

        rows.append({
            "gene": name, "chrom": chrom, "start": start, "end": end,
            "mean_depth": round(depth, 2) if depth == depth else None,
            "n_sites": len(sites.get(window, [])),
            "n_usable": n_het + n_hom,
            "n_het": n_het, "n_hom": n_hom,
            "n_lowdepth": n_lowdepth, "n_other_allele": n_other,
            "expected_het_fraction": round(expected_het_fraction, 3),
            "observed_het_fraction": (round(n_het / (n_het + n_hom), 3)
                                      if (n_het + n_hom) else None),
            "median_dev_from_half": (round(statistics.median(devs), 3) if devs else None),
            "fraction_in_central_band": round(in_band, 3) if n_het else None,
            "allelic_state": state,
            "state_statistic": stat,
        })

    bam.close()

    with open(args.out_prefix + ".genebaf.tsv", "w", encoding="utf-8") as handle:
        cols = ["gene", "chrom", "start", "end", "mean_depth", "n_sites", "n_usable",
                "n_het", "n_hom", "n_lowdepth", "n_other_allele",
                "expected_het_fraction", "observed_het_fraction",
                "median_dev_from_half", "fraction_in_central_band",
                "allelic_state", "state_statistic"]
        handle.write("\t".join(cols) + "\n")
        for row in rows:
            handle.write("\t".join("" if row[c] is None else str(row[c]) for c in cols) + "\n")

    with open(args.out_prefix + ".genebaf.sites.tsv", "w", encoding="utf-8") as handle:
        handle.write("gene\tchrom\tpos\tref_reads\talt_reads\talt_fraction\n")
        for r in per_site_rows:
            handle.write("\t".join(str(x) for x in r) + "\n")

    summary = {
        "sample": args.sample,
        "n_windows": len(windows),
        "parameters": {k: getattr(args, k) for k in
                       ("min_depth", "min_alt", "min_mapq", "min_baseq",
                        "min_hets", "loh_ratio", "min_af", "max_af", "af_field")},
        "states": {s: sum(1 for r in rows if r["allelic_state"] == s)
                   for s in ("balanced", "imbalanced", "loh", "insufficient")},
        "loh_genes": [r["gene"] for r in rows if r["allelic_state"] == "loh"],
        "imbalanced_genes": [r["gene"] for r in rows if r["allelic_state"] == "imbalanced"],
    }
    with open(args.out_prefix + ".genebaf.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"sample      : {args.sample}")
    print(f"windows     : {len(windows)}")
    print(f"states      : {summary['states']}")
    print(f"LOH         : {', '.join(summary['loh_genes']) or 'none'}")
    print(f"imbalanced  : {', '.join(summary['imbalanced_genes']) or 'none'}")
    print(f"written     : {args.out_prefix}.genebaf.tsv / .sites.tsv / .json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
