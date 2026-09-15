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
import math
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


def window_depth(bam, chrom, start, end, bin_size=0):
    """Mean depth across a window, and optionally that depth in fixed bins.

    count_coverage already builds a per-base count across the whole window;
    the original version summed it and threw the profile away. Binning the
    same arrays costs no extra read of the BAM and gives a per-gene depth
    track at a resolution the chromosome-level segment table cannot reach:
    CDKN2A is 128 kb, a fraction of one ichorCNA bin, and sits at log2 -0.90
    inside a 9p the segment caller reports as neutral.

    Returns (mean_depth, bins), where bins is a list of
    (bin_start, bin_end, mean_depth) in reference coordinates and is empty
    when bin_size is 0. The last bin of a window is short whenever the window
    length is not a multiple of bin_size; its real span is given by its own
    coordinates, so no caller needs to assume a uniform width.

    No mapping-quality filter is applied here. That matches count_coverage's
    default read_callback, which drops unmapped, secondary, QC-fail and
    duplicate reads but does not consult MAPQ, and it is the behaviour the
    reported mean has always had -- so the bins average back to exactly the
    mean this function returned before.
    """
    length = end - start
    try:
        cov = bam.count_coverage(chrom, start, end, quality_threshold=0)
    except ValueError:
        return float("nan"), []
    if not length:
        return float("nan"), []
    if not bin_size:
        return sum(sum(c) for c in cov) / length, []

    bins = []
    total = 0
    for offset in range(0, length, bin_size):
        stop = min(offset + bin_size, length)
        counts = sum(sum(c[offset:stop]) for c in cov)
        total += counts
        bins.append((start + offset, start + stop, counts / (stop - offset)))
    return total / length, bins


def classify(n_het, n_hom_usable, expected_het_fraction, devs, in_band,
             min_hets, loh_ratio):
    """Allelic state of one window from its heterozygous allele fractions.

    This does NOT call LOH. A depleted heterozygous count inside a single gene
    body is far more often homozygosity across one common haplotype than loss
    of heterozygosity: the sites in a 300 kb window sit in a handful of linkage
    blocks and are not independent observations. That test is done per arm
    instead, where the windows pooled do sit in different blocks.

    What is decided here is whether the heterozygous sites that exist are
    balanced. That reads the allele fractions directly and does not depend on
    how many hets there are.
    """
    if n_het < min_hets:
        return "insufficient", None
    median_dev = statistics.median(devs)
    # Balanced hets scatter around 0.5 with a spread set by depth; at 20x the
    # binomial SD is about 0.11, so a median deviation past 0.15 combined with
    # most hets outside the central band is a split distribution, not noise.
    if median_dev >= 0.15 and in_band < 0.35:
        return "imbalanced", round(median_dev, 3)
    return "balanced", round(median_dev, 3)


def arm_of(chrom, start, end, arm_bounds):
    """Arm label for a window, or None when it straddles the centromere."""
    if chrom not in arm_bounds:
        return None
    chrom_start, cen_start, cen_end, chrom_end = arm_bounds[chrom]
    mid = (start + end) // 2
    if mid < cen_start:
        return chrom[3:] + "p"
    if mid >= cen_end:
        return chrom[3:] + "q"
    return None


def read_arm_bounds(path):
    """{chrom: (start, cen_start, cen_end, end)} from a cytoband file."""
    spans, acen = {}, {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = canonical(parts[0])
            if chrom is None:
                continue
            start, end, stain = int(parts[1]), int(parts[2]), parts[4]
            lo, hi = spans.get(chrom, (start, end))
            spans[chrom] = (min(lo, start), max(hi, end))
            if stain == "acen":
                a, b = acen.get(chrom, (start, end))
                acen[chrom] = (min(a, start), max(b, end))
    bounds = {}
    for chrom, (lo, hi) in spans.items():
        if chrom in acen:
            bounds[chrom] = (lo, acen[chrom][0], acen[chrom][1], hi)
    return bounds


def binomial_tail(k, n, p):
    """P(X <= k) for X ~ Binomial(n, p). Exact, no scipy."""
    if n <= 0:
        return 1.0
    p = min(max(p, 1e-12), 1 - 1e-12)
    log_p, log_q = math.log(p), math.log(1 - p)
    total = 0.0
    log_coef = 0.0  # log C(n, 0)
    for i in range(0, min(k, n) + 1):
        if i > 0:
            log_coef += math.log((n - i + 1) / i)
        total += math.exp(log_coef + i * log_p + (n - i) * log_q)
    return min(1.0, total)


def arm_loh(rows, arm_bounds, min_windows, min_usable, max_p):
    """Pool windows per arm and test the heterozygous count against expectation.

    Returns {arm: {...}}. An arm covered by fewer than min_windows windows gets
    no verdict: the whole reason the per-gene test failed is that one window is
    one or two independent observations, and repeating that at arm scale would
    repeat the mistake.
    """
    pooled = {}
    for row in rows:
        arm = arm_of(row["chrom"], row["start"], row["end"], arm_bounds)
        if arm is None or row["n_usable"] == 0:
            continue
        entry = pooled.setdefault(arm, {"chrom": row["chrom"], "n_windows": 0,
                                        "n_het": 0, "n_usable": 0,
                                        "expected": [], "genes": []})
        entry["n_windows"] += 1
        entry["n_het"] += row["n_het"]
        entry["n_usable"] += row["n_usable"]
        entry["expected"].append(row["expected_het_fraction"] * row["n_usable"])
        entry["genes"].append(row["gene"])

    out = {}
    for arm, entry in sorted(pooled.items()):
        n, k = entry["n_usable"], entry["n_het"]
        expected_fraction = (sum(entry["expected"]) / n) if n else 0.0
        observed_fraction = (k / n) if n else 0.0
        ratio = (observed_fraction / expected_fraction) if expected_fraction else None

        if entry["n_windows"] < min_windows:
            verdict, pval = "not_assessed_one_window", None
        elif n < min_usable:
            verdict, pval = "not_assessed_too_few_sites", None
        else:
            pval = binomial_tail(k, n, expected_fraction)
            verdict = "loh" if pval < max_p and (ratio is not None and ratio < 0.5) \
                else "retained"

        out[arm] = {
            "chrom": entry["chrom"],
            "n_windows": entry["n_windows"],
            "genes": entry["genes"],
            "n_usable": n,
            "n_het": k,
            "expected_het_fraction": round(expected_fraction, 3),
            "observed_het_fraction": round(observed_fraction, 3),
            "ratio": round(ratio, 3) if ratio is not None else None,
            "p_value": (None if pval is None else float(f"{pval:.3g}")),
            "verdict": verdict,
        }
    return out


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
    parser.add_argument("--bin-size", type=int, default=1000,
                        help="width in bases of the per-gene depth bins. 0 "
                             "skips the bins file entirely [1000]")
    parser.add_argument("--min-hets", type=int, default=8,
                        help="hets required for a balanced/imbalanced call [8]")
    parser.add_argument("--cytobands",
                        help="cytoband file for the alignment reference. Without "
                             "it no arm-level LOH is reported.")
    parser.add_argument("--min-arm-windows", type=int, default=2,
                        help="panel windows an arm needs before LOH is tested. "
                             "One window is one or two linkage blocks, which is "
                             "what made the per-gene test unusable [2]")
    parser.add_argument("--min-arm-sites", type=int, default=100,
                        help="usable sites an arm needs before LOH is tested [100]")
    parser.add_argument("--max-arm-p", type=float, default=0.001,
                        help="binomial tail probability below which an arm is "
                             "called LOH [0.001]")
    parser.add_argument("--loh-ratio", type=float, default=0.25,
                        help="retained for compatibility; no longer used for "
                             "per-window calls")
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
    bin_rows = []

    for window in windows:
        chrom, start, end, name = window
        depth, depth_bins = window_depth(bam, chrom, start, end, args.bin_size)
        for bin_start, bin_end, bin_depth in depth_bins:
            bin_rows.append((name, chrom, bin_start, bin_end, bin_depth))
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

    arms = {}
    if args.cytobands:
        arms = arm_loh(rows, read_arm_bounds(args.cytobands),
                       args.min_arm_windows, args.min_arm_sites, args.max_arm_p)
        with open(args.out_prefix + ".armloh.tsv", "w", encoding="utf-8") as handle:
            handle.write("arm\tchrom\tn_windows\tn_usable\tn_het\t"
                         "expected_het_fraction\tobserved_het_fraction\tratio\t"
                         "p_value\tverdict\tgenes\n")
            for arm, v in arms.items():
                handle.write("\t".join(str(x) for x in (
                    arm, v["chrom"], v["n_windows"], v["n_usable"], v["n_het"],
                    v["expected_het_fraction"], v["observed_het_fraction"],
                    "" if v["ratio"] is None else v["ratio"],
                    "" if v["p_value"] is None else v["p_value"],
                    v["verdict"], ",".join(v["genes"]))) + "\n")

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

    # Depth per bin, with a log2 ratio against the median bin across the whole
    # panel for this sample. The median is over bins and not over windows, so
    # that one 5 Mb window cannot set the baseline by itself; it is recorded in
    # the JSON so a reader can renormalise without rerunning. A bin at zero
    # depth gets an empty log2 rather than negative infinity.
    normaliser = statistics.median(b[4] for b in bin_rows) if bin_rows else 0.0
    if bin_rows:
        with open(args.out_prefix + ".genebaf.bins.tsv", "w", encoding="utf-8") as handle:
            handle.write("gene\tchrom\tstart\tend\tmean_depth\tlog2_ratio\n")
            for name, chrom, bin_start, bin_end, bin_depth in bin_rows:
                if normaliser > 0 and bin_depth > 0:
                    log2 = format(math.log2(bin_depth / normaliser), ".3f")
                else:
                    log2 = ""
                handle.write(f"{name}\t{chrom}\t{bin_start}\t{bin_end}\t"
                             f"{bin_depth:.2f}\t{log2}\n")

    summary = {
        "sample": args.sample,
        "n_windows": len(windows),
        "parameters": {k: getattr(args, k) for k in
                       ("min_depth", "min_alt", "min_mapq", "min_baseq",
                        "min_hets", "loh_ratio", "min_af", "max_af", "af_field")},
        "states": {s: sum(1 for r in rows if r["allelic_state"] == s)
                   for s in ("balanced", "imbalanced", "insufficient")},
        "imbalanced_genes": [r["gene"] for r in rows if r["allelic_state"] == "imbalanced"],
        "arms": arms,
        "loh_arms": [a for a, v in arms.items() if v["verdict"] == "loh"],
        "bin_size": args.bin_size,
        "n_bins": len(bin_rows),
        "depth_normaliser": round(normaliser, 2) if bin_rows else None,
    }
    with open(args.out_prefix + ".genebaf.json", "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    print(f"sample      : {args.sample}")
    print(f"windows     : {len(windows)}")
    print(f"states      : {summary['states']}")
    if args.cytobands:
        assessed = sum(1 for v in arms.values() if v["verdict"] in ("loh", "retained"))
        print(f"arms tested : {assessed} of {len(arms)} covered")
        print(f"LOH arms    : {', '.join(summary['loh_arms']) or 'none'}")
    else:
        print("LOH arms    : not assessed (no --cytobands)")
    print(f"imbalanced  : {', '.join(summary['imbalanced_genes']) or 'none'}")
    if bin_rows:
        print(f"depth bins  : {len(bin_rows)} of {args.bin_size} bp, "
              f"normalised to {normaliser:.2f}x")
        print(f"written     : {args.out_prefix}.genebaf.tsv / .sites.tsv / "
              f".bins.tsv / .json")
    else:
        print(f"written     : {args.out_prefix}.genebaf.tsv / .sites.tsv / .json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
