#!/usr/bin/env python3
"""
mmbaf.py

Per-segment allelic imbalance from the off-target reads of an ONT adaptive
sampling run, joined to the segments produced by mmkaryo.py.

Why this exists
---------------
Read depth alone cannot separate tumour purity from tumour ploidy. Doubling
ploidy and adjusting purity maps every copy state onto another one, so a sample
fitting (purity 0.30, ploidy 2.0) usually fits (purity 0.43, ploidy 4.0) just as
well. mmkaryo.py handles that with two heuristics -- prefer the lower ploidy,
and anchor the baseline at the mode of the depth distribution -- and those
heuristics fix opposite cases and cannot both be right. In a heavily
hyperdiploid genome the modal peak is itself trisomic, and the whole karyotype
inverts into a genome of apparent losses.

Allele fractions settle it outright. A disomic region sits at a minor allele
fraction of 0.50 whatever the purity. A trisomic region at 95% purity sits at
0.34. Copy-neutral LOH sits near 0 while its depth stays flat, which depth can
never see at all.

Method
------
Per-site genotyping is hopeless at roughly 2.5x. Aggregation over a segment is
not. For every common biallelic SNP in the segment the reference and alternate
read counts are tallied, and a single minor allele fraction theta is fitted for
the whole segment by maximum likelihood, marginalising each site over
homozygous-reference, heterozygous and homozygous-alternate status using its
population allele frequency as the prior. Thousands of weakly informative sites
combine into one well-determined number.

The population frequency should match the cohort. The bundled 1000 Genomes VCF
carries per-population fields, so --af-field SAS_AF is the correct prior for a
South Asian cohort and is the default here; it sharpens theta relative to the
global AF.

Copy number is then refitted jointly against depth and allele fraction, which
removes the ambiguity rather than guessing past it.

Outputs
-------
  <prefix>.segments.baf.tsv   segments with theta, CI, site counts, allelic state
  <prefix>.baf_windows.tsv    windowed theta, for plotting as points
  <prefix>.baf.json           refitted purity/ploidy and the baseline verdict

Dependencies: numpy, pysam. No pandas, no scipy.
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
from collections import OrderedDict, defaultdict

import numpy
import pysam

from multiprocessing import Pool


CANONICAL_CHROMS = ["chr%d" % (i + 1) for i in range(22)] + ["chrX", "chrY"]

T2T_ACCESSIONS = [
    "NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1",
    "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1",
    "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1",
    "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1",
    "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1",
]

CHROM_ALIASES = {}
for _i, _canonical in enumerate(CANONICAL_CHROMS):
    CHROM_ALIASES[_canonical] = _canonical
    CHROM_ALIASES[_canonical[3:]] = _canonical
    CHROM_ALIASES[T2T_ACCESSIONS[_i]] = _canonical

BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


def canonical_chrom(name):
    return CHROM_ALIASES.get(name)


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def load_segments(path):
    """Read the segments.tsv emitted by mmkaryo.py."""
    segments = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        required = ["chromosome", "start", "end", "log2_ratio"]
        for name in required:
            if name not in index:
                raise SystemExit("ERROR: %s lacks a '%s' column. Is this an "
                                 "mmkaryo segments.tsv?" % (path, name))
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            chrom = canonical_chrom(fields[index["chromosome"]])
            if chrom is None:
                continue

            def get(name, cast=float, default=float("nan")):
                if name not in index:
                    return default
                raw = fields[index[name]]
                if raw in ("", "NA"):
                    return default
                try:
                    return cast(raw)
                except ValueError:
                    return default

            segments.append({
                "chrom": chrom,
                "start": int(fields[index["start"]]),
                "end": int(fields[index["end"]]),
                "arm": fields[index["arm"]] if "arm" in index else "",
                "n_bins": int(get("n_bins", float, 0)),
                "log2": float(fields[index["log2_ratio"]]),
                "ci_low": get("ci_low"),
                "ci_high": get("ci_high"),
                "depth_cn": get("copy_number", float, float("nan")),
            })
    return segments


def load_sites(path, af_field, fallback_field="AF", min_af=0.0, max_af=1.0):
    """Read common biallelic SNPs. Returns {chrom: (pos, ref_idx, alt_idx, af)}.

    Positions are converted to 0-based to match pysam coordinates.
    """
    per_chrom = defaultdict(lambda: ([], [], [], []))
    opener = gzip.open if path.endswith(".gz") else open
    used_fallback = 0
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t", 8)
            if len(fields) < 8:
                continue
            chrom = canonical_chrom(fields[0])
            if chrom is None:
                continue
            ref, alt = fields[3], fields[4]
            if len(ref) != 1 or len(alt) != 1:
                continue
            ref_index = BASE_INDEX.get(ref.upper())
            alt_index = BASE_INDEX.get(alt.upper())
            if ref_index is None or alt_index is None:
                continue

            info = fields[7]
            af = None
            for field_name in (af_field, fallback_field):
                if field_name is None:
                    continue
                token = field_name + "="
                start = info.find(";" + token)
                if start >= 0:
                    start += 1
                elif info.startswith(token):
                    start = 0
                else:
                    continue
                end = info.find(";", start)
                text = info[start + len(token):] if end < 0 else info[start + len(token):end]
                try:
                    af = float(text.split(",")[0])
                except ValueError:
                    af = None
                if af is not None:
                    if field_name != af_field:
                        used_fallback += 1
                    break
            if af is None or not (min_af <= af <= max_af):
                continue

            positions, refs, alts, afs = per_chrom[chrom]
            positions.append(int(fields[1]) - 1)
            refs.append(ref_index)
            alts.append(alt_index)
            afs.append(af)

    packed = OrderedDict()
    for chrom in CANONICAL_CHROMS:
        if chrom not in per_chrom:
            continue
        positions, refs, alts, afs = per_chrom[chrom]
        order = numpy.argsort(numpy.array(positions))
        packed[chrom] = (
            numpy.array(positions, dtype=numpy.int64)[order],
            numpy.array(refs, dtype=numpy.int8)[order],
            numpy.array(alts, dtype=numpy.int8)[order],
            numpy.array(afs, dtype=numpy.float64)[order],
        )
    return packed, used_fallback


# ---------------------------------------------------------------------------
# Pileup
# ---------------------------------------------------------------------------

def _tally_worker(job):
    """Tally one chromosome in its own process, with its own BAM handle.

    Profiling showed pysam.count_coverage accounts for essentially all of the
    runtime -- the array handling around it was 0.5 percent. Chromosomes are
    independent, so the only meaningful speed-up is to run them concurrently.
    """
    (bam_path, chrom, contig, length, positions, ref_index, alt_index,
     chunk_size, min_baseq, min_mapq, strict_filter) = job
    bam = pysam.AlignmentFile(bam_path, "rb")
    try:
        result = tally_chromosome(
            bam, chrom, contig, length,
            (positions, ref_index, alt_index, None),
            chunk_size, min_baseq, min_mapq, False, strict_filter)
    finally:
        bam.close()
    return chrom, result


def _coverage_views(coverage):
    """Zero-copy numpy views over the four arrays pysam.count_coverage returns.

    count_coverage hands back array.array objects. Converting them with
    numpy.asarray copies the whole window -- four arrays of five million
    elements per chunk, which is far more memory traffic than the handful of
    SNP positions actually needed. frombuffer gives a view instead.
    """
    views = []
    for column in coverage:
        try:
            views.append(numpy.frombuffer(
                memoryview(column), dtype=numpy.dtype(column.typecode)))
        except (TypeError, ValueError, AttributeError):
            views.append(numpy.asarray(column))
    return views


def tally_chromosome(bam, chrom, contig_name, length, sites, chunk_size,
                     min_baseq, min_mapq, verbose=False, strict_filter=False):
    """Reference and alternate counts at every site on one chromosome.

    pysam.count_coverage is used rather than a pileup iterator: it runs in C and
    returns per-base A/C/G/T counts for a window, so the cost is one pass over
    the reads rather than a Python-level loop over 1.4 million positions.

    Two things matter for speed here. Only the SNP offsets are pulled out of
    each window rather than the whole window being copied into numpy. And the
    read filter is left to pysam's C-level 'all' callback by default; passing a
    Python callable instead makes pysam call back into the interpreter once per
    read, which on a 16 million read BAM dominates the runtime. 'all' drops
    unmapped, secondary, QC-fail and duplicate reads but not supplementary
    alignments, and applies no mapping-quality threshold -- use --strict-filter
    to reinstate those at roughly four times the cost.
    """
    positions, ref_index, alt_index, _af = sites
    n_sites = positions.shape[0]
    n_ref = numpy.zeros(n_sites, dtype=numpy.int32)
    n_alt = numpy.zeros(n_sites, dtype=numpy.int32)
    n_other = numpy.zeros(n_sites, dtype=numpy.int32)

    if strict_filter:
        callback = lambda a: not (a.is_unmapped or a.is_secondary
                                  or a.is_supplementary or a.is_duplicate
                                  or a.mapping_quality < min_mapq)
    else:
        callback = "all"

    for chunk_start in range(0, length, chunk_size):
        chunk_end = min(chunk_start + chunk_size, length)
        lo = int(numpy.searchsorted(positions, chunk_start, side="left"))
        hi = int(numpy.searchsorted(positions, chunk_end, side="left"))
        if hi <= lo:
            continue
        try:
            coverage = bam.count_coverage(
                contig_name, chunk_start, chunk_end,
                quality_threshold=min_baseq, read_callback=callback)
        except ValueError:
            continue

        views = _coverage_views(coverage)
        offsets = (positions[lo:hi] - chunk_start).astype(numpy.int64)
        # Four small gathers of n_sites elements, not four copies of the window.
        gathered = numpy.empty((4, offsets.shape[0]), dtype=numpy.int64)
        for base in range(4):
            gathered[base] = views[base][offsets]

        columns = numpy.arange(offsets.shape[0])
        totals = gathered.sum(axis=0)
        ref_here = gathered[ref_index[lo:hi].astype(numpy.int64), columns]
        alt_here = gathered[alt_index[lo:hi].astype(numpy.int64), columns]
        n_ref[lo:hi] = ref_here
        n_alt[lo:hi] = alt_here
        n_other[lo:hi] = totals - ref_here - alt_here
        if verbose:
            sys.stderr.write("    %s:%d-%d  %d sites\n"
                             % (chrom, chunk_start, chunk_end, hi - lo))
    return n_ref, n_alt, n_other


# ---------------------------------------------------------------------------
# Minor allele fraction
# ---------------------------------------------------------------------------

THETA_GRID = numpy.linspace(0.0, 0.5, 101)


def _site_loglik(k, n, p, theta, error_rate):
    """log P(k alt reads | n reads, population AF p, minor fraction theta).

    Each site is marginalised over its three possible genotypes. For a
    heterozygous site the alternate allele may sit on either parental copy, so
    the heterozygous term is a symmetric two-component mixture. The binomial
    coefficient is common to all three terms and cancels in the argmax, so it is
    dropped.
    """
    eps = error_rate
    m = n - k

    hom_ref = (eps ** k) * ((1.0 - eps) ** m)
    hom_alt = ((1.0 - eps) ** k) * (eps ** m)

    a = theta * (1.0 - eps) + (1.0 - theta) * eps
    b = 1.0 - a
    het = 0.5 * ((a ** k) * (b ** m) + (b ** k) * (a ** m))

    total = ((1.0 - p) ** 2) * hom_ref + 2.0 * p * (1.0 - p) * het + (p ** 2) * hom_alt
    return numpy.log(numpy.maximum(total, 1e-300))


def fit_theta(k, n, p, error_rate):
    """Maximum-likelihood minor allele fraction for one segment.

    Returns (theta, ci_low, ci_high, loglik_at_theta, loglik_at_half). The last
    two give a likelihood-ratio test against the balanced 0.5 hypothesis, which
    is the question that actually matters: is this segment allelically balanced
    or not.
    """
    if k.shape[0] == 0:
        nan = float("nan")
        return nan, nan, nan, nan, nan

    profile = numpy.empty(THETA_GRID.shape[0])
    for i, theta in enumerate(THETA_GRID):
        profile[i] = float(numpy.sum(_site_loglik(k, n, p, theta, error_rate)))

    best = int(numpy.argmax(profile))
    theta_hat = float(THETA_GRID[best])
    peak = float(profile[best])

    # 95 percent profile-likelihood interval: where the log-likelihood falls
    # 1.92 below its maximum.
    inside = numpy.nonzero(profile >= peak - 1.92)[0]
    ci_low = float(THETA_GRID[inside[0]]) if inside.size else float("nan")
    ci_high = float(THETA_GRID[inside[-1]]) if inside.size else float("nan")

    balanced = float(numpy.sum(_site_loglik(k, n, p, 0.5, error_rate)))
    return theta_hat, ci_low, ci_high, peak, balanced


# ---------------------------------------------------------------------------
# Joint copy number model
# ---------------------------------------------------------------------------

def weighted_median(values, weights):
    """Length-weighted median."""
    order = numpy.argsort(values)
    values, weights = values[order], weights[order]
    cumulative = numpy.cumsum(weights)
    if cumulative[-1] <= 0:
        return float("nan")
    return float(values[int(numpy.searchsorted(cumulative, cumulative[-1] / 2.0))])


def allelic_states(max_cn):
    """All (total, major, minor) integer states up to max_cn."""
    states = []
    for total in range(0, max_cn + 1):
        for minor in range(0, total // 2 + 1):
            states.append((total, total - minor, minor))
    return states


def predicted_log2(total_cn, purity, ploidy, reference_copies=2):
    numerator = purity * total_cn + (1.0 - purity) * float(reference_copies)
    denominator = purity * ploidy + (1.0 - purity) * 2.0
    return math.log2(max(numerator, 1e-9) / denominator)


def predicted_theta(total_cn, minor_cn, purity, reference_copies=2):
    """Minor allele fraction for an allelic state at a given purity.

    Normal cells contribute one copy of each allele, so they always pull theta
    back toward 0.5. That is what makes theta purity-dependent for imbalanced
    states while leaving balanced states pinned at exactly 0.5.
    """
    normal_minor = reference_copies / 2.0
    numerator = purity * minor_cn + (1.0 - purity) * normal_minor
    denominator = purity * total_cn + (1.0 - purity) * float(reference_copies)
    if denominator <= 0:
        return float("nan")
    return min(numerator / denominator, 0.5)


def joint_fit(segments, max_cn=8, fixed_purity=None,
              purity_steps=96, ploidy_steps=121, theta_weight=1.0,
              log2_tolerance=0.05, theta_tolerance=0.03,
              tolerance=0.02, parsimony=True):
    """Refit purity and ploidy against depth and allele fraction together.

    Segments are scored against fixed tolerances, NOT against their own
    standard errors. That distinction turned out to matter enormously. A 67 Mb
    segment carries 667 bins and 12,000 informative SNPs, so its statistical
    error on log2 is about 0.003 and on theta about 0.009 -- far smaller than
    the residual any single-clone integer model leaves on real data, where
    subclonality means segments genuinely sit between states. Scoring in SE
    units therefore made the chi-square dominated by model misfit on the
    largest segments, and the optimiser chased it to nonsense: purity 0.23 on a
    sample whose depth-only fit was a well-identified 0.95.

    Fixed tolerances keep every segment's contribution bounded and make the
    score interpretable, which is the same reasoning that led mmkaryo.py to
    score in units of the copy-state gap rather than in standard errors.
    Segment length still sets the weight.
    """
    states = allelic_states(max_cn)
    purities = (numpy.array([min(max(fixed_purity, 0.01), 0.99)])
                if fixed_purity is not None
                else numpy.linspace(0.05, 0.99, purity_steps))
    ploidies = numpy.linspace(1.5, 4.5, ploidy_steps)

    usable = [s for s in segments
              if math.isfinite(s["log2"]) and s["chrom"] not in ("chrX", "chrY")]
    if not usable:
        raise SystemExit("ERROR: no autosomal segments with finite log2.")

    weights = numpy.array([max(s["n_bins"], 1) for s in usable], dtype=float)
    weights /= weights.sum()

    log2_values = numpy.array([s["log2"] for s in usable])
    theta_values = numpy.array([s.get("theta", float("nan")) for s in usable])
    has_theta = numpy.isfinite(theta_values)

    results = []
    for purity in purities:
        for ploidy in ploidies:
            expected_l = numpy.array([predicted_log2(t, purity, ploidy)
                                      for t, _a, _b in states])
            expected_t = numpy.array([predicted_theta(t, b, purity)
                                      for t, _a, b in states])

            z_l = (numpy.abs(log2_values[:, None] - expected_l[None, :])
                   / log2_tolerance)
            distance = z_l ** 2
            z_t = (numpy.abs(theta_values[:, None] - expected_t[None, :])
                   / theta_tolerance)
            distance = numpy.where(
                has_theta[:, None],
                distance + theta_weight * numpy.nan_to_num(z_t, nan=0.0) ** 2,
                distance)

            per_segment = numpy.sqrt(numpy.min(distance, axis=1))
            # Cap so that a handful of badly-fitting segments, which every real
            # tumour has, cannot drag the whole fit.
            numpy.clip(per_segment, 0.0, 3.0, out=per_segment)
            results.append((float(numpy.sum(weights * per_segment)),
                            float(purity), float(ploidy)))

    results.sort(key=lambda r: r[0])
    global_score, global_purity, global_ploidy = results[0]

    # Doubling every allelic state preserves both the depth ratio and the
    # allele fraction exactly, so whole-genome doubling is invisible to this
    # data and to any other. ASCAT and ABSOLUTE break it with a ploidy prior
    # rather than from evidence; mmkaryo.py does the same, and omitting it here
    # sent a purity-pinned fit to ploidy 3.95 on segments built at 2.0.
    if parsimony:
        cutoff = global_score + tolerance
        candidates = [r for r in results
                      if math.isfinite(r[0]) and r[0] <= cutoff]
        candidates.sort(key=lambda r: (r[2], r[0]))
        score, purity, ploidy = candidates[0] if candidates else results[0]
    else:
        score, purity, ploidy = results[0]

    return {"score": score, "purity": purity, "ploidy": ploidy,
            "global_ploidy": global_ploidy, "global_score": global_score,
            "ambiguous": abs(global_ploidy - ploidy) > 0.25,
            "fixed_purity": fixed_purity is not None,
            "n_segments_with_theta": int(numpy.count_nonzero(has_theta)),
            "n_segments": len(usable)}


def assign_state(log2_value, theta, purity, ploidy, max_cn=8,
                 reference_copies=2):
    """Best (total, major, minor) state for one segment under a fitted model."""
    best = None
    for total, major, minor in allelic_states(max_cn):
        d = abs(log2_value - predicted_log2(total, purity, ploidy,
                                            reference_copies))
        cost = (d / 0.05) ** 2
        if math.isfinite(theta):
            t = predicted_theta(total, minor, purity, reference_copies)
            cost += ((theta - t) / 0.03) ** 2
        if best is None or cost < best[0]:
            best = (cost, total, major, minor)
    return best[1], best[2], best[3]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_fai(path):
    lengths, contigs = OrderedDict(), {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            chrom = canonical_chrom(fields[0])
            if chrom is None:
                continue
            lengths[chrom] = int(fields[1])
            contigs[chrom] = fields[0]
    return lengths, contigs


def main():
    parser = argparse.ArgumentParser(
        description="Per-segment allelic imbalance from off-target adaptive "
                    "sampling reads")
    parser.add_argument("bam", help="Coordinate-sorted, indexed BAM on T2T-CHM13v2.0")
    parser.add_argument("fai", help="FASTA index for the alignment reference")
    parser.add_argument("vcf", help="Common biallelic SNP VCF in the same "
                                    "coordinates (e.g. SAVANA's "
                                    "t2t_g1000_biallelic_AF0.25-0.75.vcf.gz)")
    parser.add_argument("segments", help="segments.tsv from mmkaryo.py")
    parser.add_argument("prefix", help="Output path prefix")
    parser.add_argument("--af-field", default="SAS_AF",
                        help="INFO field to use as the population allele "
                             "frequency prior. Match it to the cohort: SAS_AF "
                             "for South Asian, EUR_AF, AFR_AF, EAS_AF, AMR_AF, "
                             "or AF for the global value (default: SAS_AF)")
    parser.add_argument("--min-depth", type=int, default=2,
                        help="Minimum reads at a site to use it (default: 2). "
                             "A single read carries no allelic information.")
    parser.add_argument("--max-depth", type=int, default=60,
                        help="Skip sites above this depth, which are usually "
                             "collapsed repeats (default: 60)")
    parser.add_argument("--min-baseq", type=int, default=10,
                        help="Minimum base quality (default: 10)")
    parser.add_argument("--min-mapq", type=int, default=1,
                        help="Minimum mapping quality (default: 1)")
    parser.add_argument("--max-other-fraction", type=float, default=0.2,
                        help="Reject a site if more than this fraction of its "
                             "reads carry neither the reference nor the "
                             "alternate base (default: 0.2)")
    parser.add_argument("--error-rate", type=float, default=0.02,
                        help="Per-base error rate for the genotype model "
                             "(default: 0.02)")
    parser.add_argument("--window", type=int, default=1000000,
                        help="Window size for the plotting-level allele "
                             "fraction track (default: 1000000)")
    parser.add_argument("--min-window-sites", type=int, default=40,
                        help="Minimum informative sites to report a window "
                             "theta (default: 40)")
    parser.add_argument("--min-sites", type=int, default=100,
                        help="Minimum informative sites to report a theta for "
                             "a segment (default: 100)")
    parser.add_argument("--purity", type=float, default=None,
                        help="Fix tumour purity instead of fitting it")
    parser.add_argument("--depth-purity", type=float, default=None,
                        help="Purity from mmkaryo's depth-only fit, used to "
                             "pin the ploidy refit when no flow purity is "
                             "available")
    parser.add_argument("--max-cn", type=int, default=8)
    parser.add_argument("--strict-filter", action="store_true",
                        help="Exclude supplementary alignments and apply "
                             "--min-mapq during the pileup. Correct but roughly "
                             "four times slower, because pysam must call back "
                             "into Python once per read. Off by default; the "
                             "site-level depth cap and --max-other-fraction "
                             "already reject the regions this would catch.")
    parser.add_argument("--chunk", type=int, default=5000000,
                        help="Pileup window size in bp (default: 5000000)")
    parser.add_argument("--jobs", type=int, default=12,
                        help="Chromosomes to tally concurrently. The pileup is "
                             "the whole runtime and chromosomes are "
                             "independent, so this is what controls wall clock "
                             "(default: 12)")
    parser.add_argument("--threads", type=int, default=2,
                        help="BAM decompression threads per worker (default: 2)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    output_dir = os.path.dirname(os.path.abspath(args.prefix))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    sys.stderr.write("Reading segments...\n")
    segments = load_segments(args.segments)
    sys.stderr.write("  %d segments\n" % len(segments))

    sys.stderr.write("Reading reference index...\n")
    lengths, contigs = load_fai(args.fai)

    sys.stderr.write("Reading SNP sites (%s prior)...\n" % args.af_field)
    sites, used_fallback = load_sites(args.vcf, args.af_field)
    total_sites = sum(v[0].shape[0] for v in sites.values())
    sys.stderr.write("  %d sites on %d contigs\n" % (total_sites, len(sites)))
    if used_fallback:
        sys.stderr.write("  note: %d sites lacked %s and fell back to AF\n"
                         % (used_fallback, args.af_field))
    if total_sites == 0:
        raise SystemExit("ERROR: no usable sites. Check the VCF and --af-field.")

    bam = pysam.AlignmentFile(args.bam, "rb", threads=args.threads)
    if not bam.has_index():
        raise SystemExit("ERROR: BAM index required. Run samtools index.")

    jobs = []
    for chrom in sites:
        if chrom not in lengths:
            continue
        contig = contigs[chrom]
        if contig not in bam.references:
            alternative = chrom if chrom in bam.references else chrom[3:]
            if alternative not in bam.references:
                continue
            contig = alternative
        positions, ref_index, alt_index, _af = sites[chrom]
        jobs.append((os.path.abspath(args.bam), chrom, contig, lengths[chrom],
                     positions, ref_index, alt_index, args.chunk,
                     args.min_baseq, args.min_mapq, args.strict_filter))
    bam.close()

    # Longest chromosomes first, so the tail of the run is not one straggler.
    jobs.sort(key=lambda j: -j[3])
    workers = max(1, min(args.jobs, len(jobs)))
    sys.stderr.write("Tallying allele counts (%d chromosomes, %d workers)...\n"
                     % (len(jobs), workers))

    tallies = {}
    if workers == 1:
        for job in jobs:
            chrom, result = _tally_worker(job)
            tallies[chrom] = result
            sys.stderr.write("  %s done\n" % chrom)
    else:
        with Pool(workers) as pool:
            for chrom, result in pool.imap_unordered(_tally_worker, jobs):
                tallies[chrom] = result
                sys.stderr.write("  %s done\n" % chrom)

    sys.stderr.write("Fitting per-segment allele fractions...\n")
    reported = 0
    for segment in segments:
        chrom = segment["chrom"]
        segment["theta"] = float("nan")
        segment["theta_ci_low"] = float("nan")
        segment["theta_ci_high"] = float("nan")
        segment["n_informative"] = 0
        segment["lrt_vs_balanced"] = float("nan")
        if chrom not in tallies:
            continue
        positions, _ref_index, _alt_index, af = sites[chrom]
        n_ref, n_alt, n_other = tallies[chrom]

        lo = int(numpy.searchsorted(positions, segment["start"], side="left"))
        hi = int(numpy.searchsorted(positions, segment["end"], side="left"))
        if hi <= lo:
            continue

        depth = n_ref[lo:hi] + n_alt[lo:hi]
        other = n_other[lo:hi]
        with numpy.errstate(divide="ignore", invalid="ignore"):
            other_fraction = other / numpy.maximum(depth + other, 1)
        keep = ((depth >= args.min_depth) & (depth <= args.max_depth)
                & (other_fraction <= args.max_other_fraction))
        segment["n_candidate"] = int(hi - lo)
        segment["n_informative"] = int(numpy.count_nonzero(keep))
        if segment["n_informative"] < args.min_sites:
            continue

        theta, ci_low, ci_high, peak, balanced = fit_theta(
            n_alt[lo:hi][keep].astype(numpy.float64),
            depth[keep].astype(numpy.float64),
            af[lo:hi][keep],
            args.error_rate)
        segment["theta"] = theta
        segment["theta_ci_low"] = ci_low
        segment["theta_ci_high"] = ci_high
        segment["lrt_vs_balanced"] = 2.0 * (peak - balanced)
        reported += 1
    sys.stderr.write("  theta estimated for %d of %d segments\n"
                     % (reported, len(segments)))

    # The joint refit is only well posed when purity is pinned. Left free it
    # can reinterpret the whole genome: (3,1) at purity 0.95 gives theta 0.339
    # and (5,1) at purity 0.33 gives 0.336, and because each segment picks its
    # own state independently there is enough freedom to match almost any
    # (purity, ploidy). On real data that produced purity 0.23 for a sample
    # whose depth-only fit was a well-identified 0.95. Per-segment theta is
    # unaffected by this and remains the trustworthy output.
    if args.purity is not None:
        sys.stderr.write("Refitting ploidy at the supplied purity...\n")
        fit = joint_fit(segments, args.max_cn, args.purity)
        sys.stderr.write("  purity %.3f (supplied), ploidy %.2f, score %.3f "
                         "(%d of %d segments contributed theta)\n"
                         % (fit["purity"], fit["ploidy"], fit["score"],
                            fit["n_segments_with_theta"], fit["n_segments"]))
    elif args.depth_purity is not None:
        sys.stderr.write("Refitting ploidy at the depth-derived purity...\n")
        fit = joint_fit(segments, args.max_cn, args.depth_purity)
        sys.stderr.write("  purity %.3f (from mmkaryo), ploidy %.2f, score "
                         "%.3f (%d of %d segments contributed theta)\n"
                         % (fit["purity"], fit["ploidy"], fit["score"],
                            fit["n_segments_with_theta"], fit["n_segments"]))
    else:
        sys.stderr.write(
            "NOT refitting purity: with purity free the joint fit is "
            "underdetermined and unreliable. Supply --purity (flow) or "
            "--depth-purity (mmkaryo's value) to pin it. Per-segment theta "
            "below is unaffected.\n")
        fit = {"purity": float("nan"), "ploidy": float("nan"),
               "score": float("nan"), "fixed_purity": False,
               "n_segments_with_theta": 0, "n_segments": len(segments),
               "not_attempted": True}

    # Is the depth baseline disomic or was it anchored on a trisomic peak?
    baseline_segments = [s for s in segments
                         if abs(s["log2"]) < 0.05 and math.isfinite(s["theta"])
                         and s["chrom"] not in ("chrX", "chrY")]
    if baseline_segments:
        # Median, not mean: copy-neutral LOH sits at log2 0 by construction, so
        # it lands in the baseline set and would drag a mean estimate down,
        # making a perfectly normal disomic baseline look imbalanced.
        weights = numpy.array([max(s["n_bins"], 1) for s in baseline_segments],
                              dtype=float)
        baseline_theta = weighted_median(
            numpy.array([s["theta"] for s in baseline_segments]), weights)
    else:
        baseline_theta = float("nan")

    if baseline_segments:
        bw = numpy.array([max(s["n_bins"], 1) for s in baseline_segments],
                         dtype=float)
        bt = numpy.array([s["theta"] for s in baseline_segments])
        total_w = bw.sum()
        baseline_spread = {
            "balanced_0.45_plus": round(float(bw[bt >= 0.45].sum() / total_w), 3),
            "imbalanced_0.28_to_0.45": round(
                float(bw[(bt >= 0.28) & (bt < 0.45)].sum() / total_w), 3),
            "strongly_imbalanced_below_0.28": round(
                float(bw[bt < 0.28].sum() / total_w), 3)}
        sys.stderr.write("  baseline theta spread: %.0f%% balanced, %.0f%% "
                         "imbalanced, %.0f%% strongly imbalanced\n"
                         % (100 * baseline_spread["balanced_0.45_plus"],
                            100 * baseline_spread["imbalanced_0.28_to_0.45"],
                            100 * baseline_spread["strongly_imbalanced_below_0.28"]))
    else:
        baseline_spread = {}

    if math.isfinite(baseline_theta):
        if baseline_theta > 0.45:
            baseline_verdict = "disomic"
        elif baseline_theta > 0.28:
            baseline_verdict = "imbalanced - baseline is probably trisomic, so "\
                               "the depth-only karyotype is shifted down by one copy"
        else:
            baseline_verdict = "strongly imbalanced - check for widespread LOH "\
                               "or a badly anchored baseline"
    else:
        baseline_verdict = "not determined"
    sys.stderr.write("  baseline theta %.3f -> %s\n"
                     % (baseline_theta, baseline_verdict))

    for segment in segments:
        if not math.isfinite(fit["purity"]):
            segment["total_cn"] = -1
            segment["major_cn"] = -1
            segment["minor_cn"] = -1
            segment["loh"] = False
            segment["cn_neutral_loh"] = False
            continue
        reference_copies = 1 if segment["chrom"] in ("chrY",) else 2
        total, major, minor = assign_state(
            segment["log2"], segment["theta"], fit["purity"], fit["ploidy"],
            args.max_cn, reference_copies)
        segment["total_cn"] = total
        segment["major_cn"] = major
        segment["minor_cn"] = minor
        segment["loh"] = bool(total > 0 and minor == 0)
        segment["cn_neutral_loh"] = bool(segment["loh"] and total == 2)

    # ---- windowed theta, for plotting -------------------------------------
    # Per-segment theta is the right statistic for calling, but a segment table
    # gives a plot nothing to scatter -- you get sparse horizontal bars instead
    # of the dot cloud that TitanCNA and PURPLE's circos show. Fixed windows
    # give the density back. At 2.5x a 1 Mb window holds roughly 145 usable
    # het sites, enough for theta to about +/-0.03, whereas a single site at
    # depth 2 or 3 carries essentially no information and would just produce
    # stripes at 0, 1/3, 1/2.
    sys.stderr.write("Fitting windowed allele fractions (%d bp windows)...\n"
                     % args.window)
    windows = []
    for chrom in tallies:
        positions, _ri, _ai, af = sites[chrom]
        n_ref, n_alt, n_other = tallies[chrom]
        depth = n_ref + n_alt
        with numpy.errstate(divide="ignore", invalid="ignore"):
            other_fraction = n_other / numpy.maximum(depth + n_other, 1)
        usable = ((depth >= args.min_depth) & (depth <= args.max_depth)
                  & (other_fraction <= args.max_other_fraction))
        if positions.shape[0] == 0:
            continue
        limit = int(positions[-1]) + args.window
        for start in range(0, limit, args.window):
            lo = int(numpy.searchsorted(positions, start, side="left"))
            hi = int(numpy.searchsorted(positions, start + args.window,
                                        side="left"))
            if hi <= lo:
                continue
            keep = usable[lo:hi]
            if int(numpy.count_nonzero(keep)) < args.min_window_sites:
                continue
            theta, ci_low, ci_high, _peak, _bal = fit_theta(
                n_alt[lo:hi][keep].astype(numpy.float64),
                depth[lo:hi][keep].astype(numpy.float64),
                af[lo:hi][keep], args.error_rate)
            windows.append({"chrom": chrom, "start": start,
                            "end": start + args.window,
                            "n_sites": int(numpy.count_nonzero(keep)),
                            "theta": theta})
    sys.stderr.write("  %d windows with theta\n" % len(windows))

    with open(args.prefix + ".baf_windows.tsv", "w") as handle:
        handle.write("chromosome\tstart\tend\tn_informative_sites\ttheta\n")
        for w in windows:
            handle.write("%s\t%d\t%d\t%d\t%s\n"
                         % (w["chrom"], w["start"], w["end"], w["n_sites"],
                            ("%.4f" % w["theta"])
                            if math.isfinite(w["theta"]) else "NA"))

    sys.stderr.write("Writing outputs...\n")
    with open(args.prefix + ".segments.baf.tsv", "w") as handle:
        handle.write("chromosome\tstart\tend\tarm\tn_bins\tsize_mb\t"
                     "log2_ratio\tn_candidate_sites\tn_informative_sites\t"
                     "theta\ttheta_ci_low\ttheta_ci_high\tlrt_vs_balanced\t"
                     "total_cn\tmajor_cn\tminor_cn\tloh\tcn_neutral_loh\n")
        for s in segments:
            def show(value, fmt="%.4f"):
                return fmt % value if math.isfinite(value) else "NA"
            handle.write("%s\t%d\t%d\t%s\t%d\t%.2f\t%.4f\t%d\t%d\t%s\t%s\t%s\t"
                         "%s\t%d\t%d\t%d\t%s\t%s\n" % (
                s["chrom"], s["start"], s["end"], s["arm"], s["n_bins"],
                (s["end"] - s["start"]) / 1e6, s["log2"],
                s.get("n_candidate", 0), s["n_informative"],
                show(s["theta"], "%.4f"), show(s["theta_ci_low"], "%.4f"),
                show(s["theta_ci_high"], "%.4f"),
                show(s["lrt_vs_balanced"], "%.1f"),
                s["total_cn"], s["major_cn"], s["minor_cn"],
                "yes" if s["loh"] else "no",
                "yes" if s["cn_neutral_loh"] else "no"))

    cn_loh = [s for s in segments if s["cn_neutral_loh"]
              and (s["end"] - s["start"]) >= 3000000]
    summary = {
        "bam": os.path.abspath(args.bam),
        "segments_input": os.path.abspath(args.segments),
        "af_field": args.af_field,
        "parameters": {
            "min_depth": args.min_depth, "max_depth": args.max_depth,
            "min_baseq": args.min_baseq, "min_mapq": args.min_mapq,
            "error_rate": args.error_rate, "min_sites": args.min_sites,
            "max_other_fraction": args.max_other_fraction,
        },
        "sites_loaded": total_sites,
        "segments_with_theta": reported,
        "fit": {"purity": round(fit["purity"], 4),
                "ploidy": round(fit["ploidy"], 3),
                "score": round(fit["score"], 4),
                "purity_supplied": fit["fixed_purity"],
                "segments_contributing_theta": fit["n_segments_with_theta"]},
        "baseline": {"theta": (round(baseline_theta, 4)
                               if math.isfinite(baseline_theta) else None),
                     "verdict": baseline_verdict,
                     "spread": baseline_spread,
                     "n_segments": len(baseline_segments)},
        "cn_neutral_loh": [
            {"region": "%s:%d-%d" % (s["chrom"], s["start"], s["end"]),
             "size_mb": round((s["end"] - s["start"]) / 1e6, 2),
             "theta": round(s["theta"], 4)} for s in cn_loh],
    }
    with open(args.prefix + ".baf.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    if cn_loh:
        sys.stderr.write("  copy-neutral LOH in %d segments >= 3 Mb\n"
                         % len(cn_loh))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
