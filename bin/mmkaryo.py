#!/usr/bin/env python3
"""
mmkaryo.py

Digital karyotyping and copy-number analysis from Oxford Nanopore adaptive
sampling BAM files, tuned for plasma cell neoplasms.

Design notes
------------
The input is an adaptive-sampling BAM in which a small fraction of the genome
(the enrichment panel) is deeply covered and the remainder is covered by the
short "rejected" reads that adaptive sampling produces when it ejects an
off-target molecule. Those rejected reads are numerous even though their base
coverage is low, so copy number is inferred from READ COUNTS PER BIN rather
than from nucleotide depth.

The enrichment panel is deliberately excluded from the analysis via the
exclusion mask, so panel boundaries cannot create artificial copy-number
changepoints. Pass the panel BED to --exclude alongside repeat and satellite
annotation.

Pipeline stages
---------------
  1. Bin the genome at a fixed bin size.
  2. Build a fine-resolution exclusion mask from one or more BED files.
  3. Stream the BAM once, counting one read per fragment at its midpoint,
     discarding reads whose midpoint falls inside a masked interval.
  4. Correct counts for the unmasked fraction of each bin, then optionally for
     GC content.
  5. Convert to log2 ratio against the autosomal median.
  6. Segment each chromosome arm by recursive binary segmentation.
  7. Fit tumour purity and ploidy by grid search over integer copy states, or
     accept an externally measured purity.
  8. Report per-segment absolute copy number, implied clonal cell fraction,
     and the detection limit given the observed noise.

Outputs
-------
  <prefix>.bins.tsv       per-bin counts and log2 ratios
  <prefix>.segments.tsv   segment table with copy number and cell fraction
  <prefix>.arms.tsv       arm-level summary (comparable to lrdk.py output)
  <prefix>.cytobands.tsv  per-cytoband summary (only with --cytobands)
  <prefix>.karyotype.json machine-readable summary including the ISCN string
  <prefix>.karyotype.png  genome-wide plot (only if matplotlib is available)

Dependencies: numpy, pysam. matplotlib is optional and only used for plotting.
No pandas, no scipy.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from collections import OrderedDict, defaultdict

import numpy
import pysam


# ---------------------------------------------------------------------------
# Reference constants (T2T-CHM13v2.0)
# ---------------------------------------------------------------------------

# RefSeq accessions for CHM13v2.0, chr1..chr22, chrX, chrY in order.
T2T_ACCESSIONS = [
    "NC_060925.1", "NC_060926.1", "NC_060927.1", "NC_060928.1", "NC_060929.1",
    "NC_060930.1", "NC_060931.1", "NC_060932.1", "NC_060933.1", "NC_060934.1",
    "NC_060935.1", "NC_060936.1", "NC_060937.1", "NC_060938.1", "NC_060939.1",
    "NC_060940.1", "NC_060941.1", "NC_060942.1", "NC_060943.1", "NC_060944.1",
    "NC_060945.1", "NC_060946.1", "NC_060947.1", "NC_060948.1",
]

CANONICAL_CHROMS = ["chr%d" % (i + 1) for i in range(22)] + ["chrX", "chrY"]

# CHM13v2.0 centromere spans, used only to split chromosomes into p and q arms.
# Acrocentric chromosomes (13, 14, 15, 21, 22) and the sex chromosomes are
# handled as single units, matching standard cytogenetic reporting.
T2T_CENTROMERES = {
    "chr1":  (121405145, 144008994),
    "chr2":  (87403406, 98345164),
    "chr3":  (90339372, 96498140),
    "chr4":  (49050159, 57435402),
    "chr5":  (46161727, 51082874),
    "chr6":  (57216843, 62665581),
    "chr7":  (55900922, 66709007),
    "chr8":  (43506737, 48466602),
    "chr9":  (40119875, 81269165),
    "chr10": (38527422, 43443011),
    "chr11": (48819609, 55216869),
    "chr12": (30997656, 38541929),
    "chr16": (32625033, 52263389),
    "chr17": (20527834, 28139193),
    "chr18": (11350790, 21125454),
    "chr19": (19876180, 30396484),
    "chr20": (25835813, 37114981),
}

# Chromosomes gained in hyperdiploid multiple myeloma.
HYPERDIPLOID_CHROMS = ["chr3", "chr5", "chr7", "chr9", "chr11",
                       "chr15", "chr19", "chr21"]

# Contig-name lookup accepting NC_ accessions, chr-prefixed names and bare
# numerals, all resolving to the canonical chr-prefixed form.
CHROM_ALIASES = {}
for _i, _canonical in enumerate(CANONICAL_CHROMS):
    CHROM_ALIASES[_canonical] = _canonical
    CHROM_ALIASES[_canonical[3:]] = _canonical
    CHROM_ALIASES[T2T_ACCESSIONS[_i]] = _canonical


def canonical_chrom(name):
    """Map any accepted contig name to chr1..chr22/chrX/chrY, else None."""
    return CHROM_ALIASES.get(name)


# ---------------------------------------------------------------------------
# Binning and masking
# ---------------------------------------------------------------------------

class Genome(object):
    """Fixed-width bins over the canonical chromosomes, with an exclusion mask.

    The mask is stored at a coarse resolution (default 1 kb blocks) rather than
    per base. At a 200 kb bin size this costs nothing in accuracy and keeps the
    whole-genome mask under a few megabytes instead of a few gigabytes.
    """

    def __init__(self, lengths, bin_size, mask_resolution):
        self.bin_size = int(bin_size)
        self.mask_resolution = int(mask_resolution)
        self.lengths = OrderedDict()
        self.n_bins = OrderedDict()
        self.mask = OrderedDict()      # True where excluded
        self.counts = OrderedDict()

        for chrom in CANONICAL_CHROMS:
            if chrom not in lengths:
                continue
            length = int(lengths[chrom])
            self.lengths[chrom] = length
            self.n_bins[chrom] = length // self.bin_size + 1
            n_blocks = length // self.mask_resolution + 1
            self.mask[chrom] = numpy.zeros(n_blocks, dtype=bool)
            self.counts[chrom] = numpy.zeros(self.n_bins[chrom], dtype=numpy.int64)

    def add_mask_bed(self, path):
        """Mark every interval in a BED file as excluded. Returns rows used."""
        used = 0
        with open(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                chrom = canonical_chrom(fields[0])
                if chrom is None or chrom not in self.mask:
                    continue
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    continue
                if end <= start:
                    continue
                lo = start // self.mask_resolution
                hi = min(end // self.mask_resolution + 1, self.mask[chrom].shape[0])
                self.mask[chrom][lo:hi] = True
                used += 1
        return used

    def mask_centromeres(self, pad=0):
        """Exclude centromere spans, which are unreliable for read counting."""
        for chrom, (start, end) in T2T_CENTROMERES.items():
            if chrom not in self.mask:
                continue
            lo = max(0, (start - pad) // self.mask_resolution)
            hi = min((end + pad) // self.mask_resolution + 1,
                     self.mask[chrom].shape[0])
            self.mask[chrom][lo:hi] = True

    def unmasked_fraction(self):
        """Fraction of each bin that survives the mask, per chromosome."""
        blocks_per_bin = self.bin_size // self.mask_resolution
        fractions = OrderedDict()
        for chrom in self.lengths:
            mask = self.mask[chrom]
            n_bins = self.n_bins[chrom]
            padded = numpy.zeros(n_bins * blocks_per_bin, dtype=bool)
            padded[:mask.shape[0]] = mask[:padded.shape[0]]
            # Blocks beyond the contig end count as masked.
            if mask.shape[0] < padded.shape[0]:
                padded[mask.shape[0]:] = True
            reshaped = padded.reshape(n_bins, blocks_per_bin)
            fractions[chrom] = 1.0 - reshaped.mean(axis=1)
        return fractions


def count_reads(bam_path, genome, min_mapq, threads, read_limit=None,
                verbose=False, allow_truncation=False):
    """Stream the BAM once and count one read per fragment at its midpoint.

    Secondary and supplementary alignments are skipped so that a split read
    contributes exactly one count. lrdk.py counts every alignment record, which
    double-counts reads spanning a structural variant breakpoint.
    """
    stats = {"total": 0, "counted": 0, "masked": 0, "unplaced": 0,
             "low_mapq": 0, "secondary_supp": 0, "unmapped": 0}

    # pysam refuses extra decompression threads when ignoring truncation, so
    # the two are mutually exclusive. Threads are the better default for a
    # completed BAM; --allow-truncated trades them for the ability to run
    # against a sequencing run that is still in progress.
    if allow_truncation:
        handle = pysam.AlignmentFile(bam_path, mode="rb", ignore_truncation=True)
    else:
        handle = pysam.AlignmentFile(bam_path, mode="rb", threads=threads)

    # Fail early on a reference mismatch rather than silently returning zeros.
    recognised = sum(1 for ref in handle.references
                     if canonical_chrom(ref) is not None)
    if recognised == 0:
        handle.close()
        raise SystemExit(
            "ERROR: no BAM contig names map to T2T-CHM13v2.0 (expected NC_0609*"
            " accessions or chr-prefixed names). Check the alignment reference.")

    mask_res = genome.mask_resolution
    bin_size = genome.bin_size

    for aln in handle.fetch(until_eof=True):
        stats["total"] += 1
        if read_limit is not None and stats["total"] > read_limit:
            break
        if aln.is_unmapped:
            stats["unmapped"] += 1
            continue
        if aln.is_secondary or aln.is_supplementary:
            stats["secondary_supp"] += 1
            continue
        if aln.mapping_quality < min_mapq:
            stats["low_mapq"] += 1
            continue
        chrom = canonical_chrom(aln.reference_name)
        if chrom is None or chrom not in genome.counts:
            stats["unplaced"] += 1
            continue
        start = aln.reference_start
        end = aln.reference_end
        if end is None:
            end = start + 1
        midpoint = (start + end) // 2
        block = midpoint // mask_res
        chrom_mask = genome.mask[chrom]
        if block >= chrom_mask.shape[0] or chrom_mask[block]:
            stats["masked"] += 1
            continue
        index = midpoint // bin_size
        if index >= genome.counts[chrom].shape[0]:
            stats["unplaced"] += 1
            continue
        genome.counts[chrom][index] += 1
        stats["counted"] += 1

        if verbose and stats["total"] % 5000000 == 0:
            sys.stderr.write("  %d reads read, %d counted\n"
                             % (stats["total"], stats["counted"]))

    handle.close()
    return stats


# ---------------------------------------------------------------------------
# Normalisation
# ---------------------------------------------------------------------------

def gc_content(reference_path, genome, chrom, n_bins):
    """GC fraction of each bin, ignoring N. Returns NaN for all-N bins."""
    fasta = reference_path  # already an open pysam.FastaFile
    values = numpy.full(n_bins, numpy.nan, dtype=numpy.float64)
    length = genome.lengths[chrom]
    # pysam requires the contig name as it appears in the FASTA index.
    contig = genome.fasta_contig[chrom]
    for i in range(n_bins):
        start = i * genome.bin_size
        end = min(start + genome.bin_size, length)
        if end <= start:
            continue
        seq = fasta.fetch(contig, start, end).upper()
        gc = seq.count("G") + seq.count("C")
        at = seq.count("A") + seq.count("T")
        if gc + at > 0:
            values[i] = float(gc) / float(gc + at)
    return values


def correct_gc(counts, gc, n_strata=50):
    """Divide out the GC-dependent component of read count.

    A median-per-GC-stratum correction is used rather than a LOESS fit. It is
    robust to outliers, has no tuning parameters beyond the stratum count, and
    needs no external dependency. Bins with unknown GC are left uncorrected.
    """
    corrected = counts.astype(numpy.float64).copy()
    usable = numpy.isfinite(gc) & (counts > 0)
    if usable.sum() < n_strata * 10:
        return corrected, False

    overall = numpy.median(counts[usable])
    edges = numpy.quantile(gc[usable], numpy.linspace(0.0, 1.0, n_strata + 1))
    edges[0] -= 1e-9
    edges[-1] += 1e-9
    stratum = numpy.digitize(gc, edges) - 1

    for s in range(n_strata):
        members = usable & (stratum == s)
        if members.sum() < 10:
            continue
        local = numpy.median(counts[members])
        if local <= 0:
            continue
        corrected[members] = counts[members] * (overall / local)
    return corrected, True


def modal_baseline(values, sigma_hint=None):
    """Most common bin depth, used as the diploid reference level.

    The median is the obvious choice and it is wrong for this disease. In
    hyperdiploid myeloma a quarter to a half of the genome carries a trisomy,
    and that aneuploid mass drags the 50th percentile up into the upper flank
    of the disomic peak. On a simulated sample with 27 percent of bins trisomic
    the median baseline sat 4 percent above the true disomic level, which turns
    every normal region into an apparent loss and pushes genuine low-level
    gains below the reporting threshold.

    The mode of the depth distribution does not have that problem: extra
    aneuploid mass forms its own peak instead of shifting the dominant one.

    Returns (baseline, method, median_baseline) so the caller can report how
    far apart the two estimates were, which is itself a useful indicator of
    how much of the genome is altered.
    """
    finite = values[numpy.isfinite(values) & (values > 0)]
    median_value = float(numpy.median(finite)) if finite.shape[0] else float("nan")
    if finite.shape[0] < 200:
        return median_value, "median", median_value

    # Work in log space so that the peak width is copy-ratio invariant.
    logged = numpy.log2(finite)
    low, high = numpy.percentile(logged, [1.0, 99.0])
    if not (high > low):
        return median_value, "median", median_value

    # Kernel density by convolution of a fine histogram. The bandwidth must
    # track the per-bin noise: too narrow and the estimate chases the integer
    # granularity of read counts rather than the copy-state peak, which is how
    # an earlier histogram version of this function silently reproduced the
    # median it was meant to replace.
    if sigma_hint and sigma_hint > 0:
        bandwidth = max(sigma_hint / 2.0, (high - low) / 200.0)
    else:
        bandwidth = (high - low) / 25.0
    bandwidth = min(bandwidth, (high - low) / 6.0)

    grid_size = 1024
    counts, edges = numpy.histogram(logged, bins=grid_size, range=(low, high))
    step = edges[1] - edges[0]

    radius = max(1, int(math.ceil(3.0 * bandwidth / step)))
    offsets = numpy.arange(-radius, radius + 1) * step
    kernel = numpy.exp(-0.5 * (offsets / bandwidth) ** 2)
    kernel /= kernel.sum()
    density = numpy.convolve(counts.astype(float), kernel, mode="same")

    peak = int(numpy.argmax(density))
    centre = 0.5 * (edges[peak] + edges[peak + 1])

    # Refine by taking the median of the points under the peak.
    near = numpy.abs(logged - centre) < bandwidth
    if numpy.count_nonzero(near) >= 50:
        centre = float(numpy.median(logged[near]))
    return float(2.0 ** centre), "mode", median_value


def robust_sd_from_diffs(values):
    """Robust per-bin SD estimated from successive differences.

    Successive differencing removes any smooth copy-number structure, so the
    scatter that remains is technical noise. The MAD makes it insensitive to
    the genuine changepoints.
    """
    finite = values[numpy.isfinite(values)]
    if finite.shape[0] < 3:
        return float("nan")
    diffs = numpy.diff(finite)
    mad = numpy.median(numpy.abs(diffs - numpy.median(diffs)))
    return float(mad * 1.4826 / math.sqrt(2.0))


# ---------------------------------------------------------------------------
# Segmentation
# ---------------------------------------------------------------------------

MAX_ARC_EVALUATIONS = 4000000


def _best_arc(values, sigma, min_bins):
    """Best contiguous interval against its complement (the CBS statistic).

    Plain binary segmentation scans single split points, which fails on a small
    focal event sitting in the middle of a long arm: at the first split the
    event is diluted across one whole side and the statistic never clears
    threshold. A 15 percent subclonal 1q21 gain on a 100 Mb arm is exactly that
    case. Circular binary segmentation instead scans intervals, so a focal
    change is compared against everything outside it and is found directly.

    Returns (t, start, end). start/end delimit the candidate interval; if it
    touches an edge of the window the caller treats it as a single split.
    """
    n = values.shape[0]
    if n < 2 * min_bins or not math.isfinite(sigma) or sigma <= 0:
        return 0.0, None, None

    cumulative = numpy.concatenate(([0.0], numpy.cumsum(values)))
    total = cumulative[-1]

    lengths = list(range(min_bins, n - min_bins + 1))
    # Cap the work on very fine bin sizes by striding the interval-length grid.
    # Boundaries are refined by the recursion, so a coarse grid costs little.
    if len(lengths) * n > MAX_ARC_EVALUATIONS:
        stride = max(1, int(len(lengths) * n / MAX_ARC_EVALUATIONS))
        lengths = lengths[::stride]

    best_t, best_start, best_end = 0.0, None, None
    for inside_length in lengths:
        outside_length = n - inside_length
        if outside_length < min_bins:
            continue
        inside_sum = cumulative[inside_length:] - cumulative[:n - inside_length + 1]
        outside_sum = total - inside_sum
        difference = inside_sum / inside_length - outside_sum / outside_length
        scale = sigma * math.sqrt(1.0 / inside_length + 1.0 / outside_length)
        t = numpy.abs(difference) / scale
        index = int(numpy.argmax(t))
        if t[index] > best_t:
            best_t = float(t[index])
            best_start = index
            best_end = index + inside_length
    return best_t, best_start, best_end


def segment_values(values, sigma, min_bins, threshold):
    """Recursive circular binary segmentation. Returns (start, end) index pairs.

    Splits are accepted on statistical significance alone. There is deliberately
    no effect-size gate here, because applying one at this stage would discard
    the diluted intermediate splits that lead to real focal events. Effect size
    is applied afterwards, in merge_segments.

    The threshold rises with window length because the statistic is a maximum
    over order n^2 candidate intervals, whose null maximum grows like
    sqrt(2 log(n^2)) = 2 sqrt(log n).
    """
    boundaries = []

    def recurse(lo, hi):
        n = hi - lo
        if n < 2 * min_bins:
            boundaries.append((lo, hi))
            return
        window = values[lo:hi]
        t, start, end = _best_arc(window, sigma, min_bins)
        if start is None:
            boundaries.append((lo, hi))
            return
        effective = max(threshold, 2.0 * math.sqrt(math.log(max(n, 3))) + 0.5)
        if t < effective:
            boundaries.append((lo, hi))
            return
        # Split into up to three pieces, dropping any zero-length edge piece.
        cuts = [0, start, end, n]
        cuts = sorted(set(c for c in cuts if 0 <= c <= n))
        pieces = [(cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]
        if len(pieces) < 2 or any((b - a) < min_bins for a, b in pieces):
            boundaries.append((lo, hi))
            return
        for a, b in pieces:
            recurse(lo + a, lo + b)

    recurse(0, values.shape[0])
    boundaries.sort()
    return boundaries


def merge_segments(values, boundaries, min_delta):
    """Collapse adjacent segments whose medians differ by less than min_delta.

    Applied repeatedly until stable, because each merge changes the median that
    the next comparison is made against.
    """
    if not boundaries:
        return boundaries
    current = [list(b) for b in boundaries]
    while True:
        merged = [current[0]]
        for lo, hi in current[1:]:
            previous = merged[-1]
            a = numpy.median(values[previous[0]:previous[1]])
            b = numpy.median(values[lo:hi])
            if abs(a - b) < min_delta:
                previous[1] = hi
            else:
                merged.append([lo, hi])
        if len(merged) == len(current):
            return [tuple(m) for m in merged]
        current = merged


# ---------------------------------------------------------------------------
# Absolute copy number
# ---------------------------------------------------------------------------

def expected_log2(copy_numbers, purity, ploidy, reference_copies=2):
    """log2 depth ratio for each integer copy number under (purity, ploidy).

    ratio(n) = (rho * n + (1 - rho) * ref) / (rho * psi + (1 - rho) * 2)

    reference_copies is the germline copy number of the region: 2 everywhere
    except the sex chromosomes of a male sample, where it is 1. It enters the
    normal-cell term of the numerator. Note that it must NOT be applied by
    shifting the observed log2 instead, which would report a normal male X as
    two copies.
    """
    numerator = purity * copy_numbers + (1.0 - purity) * float(reference_copies)
    denominator = purity * ploidy + (1.0 - purity) * 2.0
    with numpy.errstate(divide="ignore"):
        return numpy.log2(numpy.maximum(numerator, 1e-9) / denominator)


def _scaled_distance(seg_log2, expected):
    """Distance to the nearest integer state, in units of the local state gap.

    Raw log2 distance must not be used as a fit score. Increasing ploidy packs
    the integer copy states closer together in log2 space, so a raw-distance
    score can always be improved by inflating ploidy until some state happens
    to sit near every observed segment. Dividing by the spacing of the two
    bracketing states removes that bias and bounds the result at 0.5, which is
    what a segment lying exactly between two states scores.
    """
    gaps = numpy.diff(expected)
    gaps = numpy.where(gaps > 1e-6, gaps, 1e-6)
    differences = numpy.abs(seg_log2[:, None] - expected[None, :])
    nearest = numpy.argmin(differences, axis=1)
    raw = differences[numpy.arange(seg_log2.shape[0]), nearest]

    left = numpy.clip(nearest - 1, 0, gaps.shape[0] - 1)
    right = numpy.clip(nearest, 0, gaps.shape[0] - 1)
    local_gap = numpy.minimum(gaps[left], gaps[right])
    return numpy.minimum(raw / local_gap, 0.5)


def fit_purity_ploidy(seg_log2, seg_weight, fixed_purity=None, max_cn=8,
                      purity_steps=192, ploidy_steps=121, tolerance=0.01,
                      parsimony=True):
    """Grid search over purity and ploidy minimising distance to integer states.

    Purity and ploidy are not jointly identifiable from read depth alone. A
    sample fitting (purity 0.30, ploidy 2.0) will usually fit (purity 0.43,
    ploidy 4.0) just as well, because doubling ploidy and adjusting purity maps
    every copy state onto another one. Allele-fraction data is what normally
    breaks the tie, and this script does not use it.

    Two things are therefore done rather than pretending the ambiguity is not
    there. Among all fits within `tolerance` of the best score, the lowest
    ploidy is preferred, on the grounds that whole-genome doubling is far less
    common in plasma cell neoplasms than hyperdiploidy. And the global best is
    reported alongside, together with the margin separating ploidy solutions,
    so an ambiguous sample is visible as such.

    Supplying --purity from post-sort flow cytometry removes the ambiguity
    outright and is the recommended mode.
    """
    copy_numbers = numpy.arange(0, max_cn + 1, dtype=numpy.float64)
    if fixed_purity is not None:
        purities = numpy.array([min(max(fixed_purity, 0.01), 0.99)])
    else:
        purities = numpy.linspace(0.05, 0.99, purity_steps)
    ploidies = numpy.linspace(1.5, 4.5, ploidy_steps)

    weight = seg_weight / float(numpy.sum(seg_weight))
    results = []
    for purity in purities:
        for ploidy in ploidies:
            expected = expected_log2(copy_numbers, purity, ploidy)
            distance = _scaled_distance(seg_log2, expected)
            results.append((float(numpy.sum(weight * distance)),
                            float(purity), float(ploidy)))
    results.sort(key=lambda r: r[0])
    global_score, global_purity, global_ploidy = results[0]

    if parsimony:
        cutoff = global_score + tolerance
        candidates = [r for r in results
                      if math.isfinite(r[0]) and r[0] <= cutoff]
        if candidates:
            candidates.sort(key=lambda r: (r[2], r[0]))
            best_score, best_purity, best_ploidy = candidates[0]
        else:
            best_score, best_purity, best_ploidy = (
                global_score, global_purity, global_ploidy)
    else:
        best_score, best_purity, best_ploidy = global_score, global_purity, global_ploidy

    margin = float("nan")
    for score, _purity, ploidy in results:
        if abs(ploidy - best_ploidy) > 0.25:
            margin = score - best_score
            break

    # Fraction of the genome that sits closer to halfway between two integer
    # states than to either of them. In a plasma cell neoplasm this is usually
    # real subclonality rather than a bad fit: no single purity value can
    # describe a sample carrying a 60 percent trisomy and a 30 percent
    # monosomy at the same time. For those segments the per-segment cell
    # fraction is the meaningful readout, not the integer copy number.
    expected = expected_log2(copy_numbers, best_purity, best_ploidy)
    residuals = _scaled_distance(seg_log2, expected)
    unexplained = float(numpy.sum(weight[residuals > 0.25]))

    return {
        "score": best_score,
        "purity": best_purity,
        "ploidy": best_ploidy,
        "margin": margin,
        "unexplained_fraction": unexplained,
        "global_score": global_score,
        "global_purity": global_purity,
        "global_ploidy": global_ploidy,
        "ambiguous": (not fixed_purity) and abs(global_ploidy - best_ploidy) > 0.25,
        "fixed_purity": fixed_purity is not None,
    }


def assign_copy_number(log2_value, purity, ploidy, max_cn=8, reference_copies=2):
    """Nearest integer copy state and its residual, given a fitted model.

    reference_copies is the expected germline copy number of the region: 2
    everywhere except the sex chromosomes of a male sample, where it is 1 so
    that a normal X is reported as one copy rather than as a deletion.
    """
    copy_numbers = numpy.arange(0, max_cn + 1, dtype=numpy.float64)
    expected = expected_log2(copy_numbers, purity, ploidy, reference_copies)
    index = int(numpy.argmin(numpy.abs(expected - log2_value)))
    return int(copy_numbers[index]), float(log2_value - expected[index])


def cell_fraction(log2_value, expected_copies=2):
    """Clonal cell fraction implied by a log2 ratio.

    Observed copies relative to the autosomal baseline are 2 * 2**log2, so a
    single-copy change in a fraction f of cells gives

        f = |2 * 2**log2 - expected_copies|

    On a diploid background (expected_copies = 2) this is the original
    formulation, 2 * (2**log2 - 1) for a gain and 2 * (1 - 2**log2) for a loss.
    Passing expected_copies = 1 makes it correct for X and Y in a male, where
    the baseline is one copy and a normal chromosome is not an event.

    These fractions are directly comparable to the percentage of nuclei scored
    positive by interphase FISH, which is what makes them useful for
    concordance work. NaN is returned when the implied fraction falls outside
    (0, 1], which is the case for changes of more than one copy.
    """
    observed_copies = 2.0 * (2.0 ** log2_value)
    value = abs(observed_copies - float(expected_copies))
    if not math.isfinite(value) or value <= 0.0 or value > 1.0:
        return float("nan")
    return float(value)


def classify_event(copy_number, log2_value, expected_copies=2):
    """Label a region from its FITTED copy number, not from the sign of log2.

    Returns (kind, fraction) where kind is "gain", "loss" or "none".

    The sign of log2 alone cannot distinguish a real event from baseline noise:
    a region at the fitted copy number still has a log2 a little either side of
    zero, and labelling that a gain or a loss produced a finding on almost
    every arm of a near-diploid genome. The integer state is the assignment the
    model actually made, so the event follows from it.

    A region sitting at its expected copy number gets no cell fraction. A cell
    fraction is the clonal fraction OF AN EVENT; reporting 0.012 for a neutral
    arm invites it to be read as a 1.2 percent subclone.
    """
    if expected_copies is None:
        # chrY in a female. There is no expected copy number because there is
        # no chromosome: the observed log2 is a mapping floor produced by
        # mismapped repetitive reads, not a measurement of anything.
        return "none", float("nan")
    if copy_number == expected_copies:
        return "none", float("nan")
    kind = "gain" if copy_number > expected_copies else "loss"
    return kind, cell_fraction(log2_value, expected_copies)


def event_expected_copies(chrom, sex_label):
    """Expected copy number for EVENT CALLING, or None where it does not apply.

    Distinct from the reference copies used by assign_copy_number, which needs
    a positive number to compute an expected log2 and therefore cannot take
    None. Here a None means the region should not be scored at all.

        autosome        2
        chrX in XY      1      a normal male X is not a deletion
        chrX in XX      2      a single X in a female IS a deletion
        chrY in XY      1
        chrY in XX      None   not scored
    """
    if chrom == "chrX":
        return 1 if sex_label == "XY" else 2
    if chrom == "chrY":
        return 1 if sex_label == "XY" else None
    return 2


def detection_limit(sigma, n_bins, n_sigma=3.0):
    """Smallest single-copy cell fraction detectable over n_bins at n_sigma."""
    if not math.isfinite(sigma) or n_bins <= 0:
        return float("nan")
    standard_error = sigma / math.sqrt(n_bins)
    return float(2.0 * (2.0 ** (n_sigma * standard_error) - 1.0))


# ---------------------------------------------------------------------------
# Arm and cytoband summaries
# ---------------------------------------------------------------------------

def infer_sex(bins, baseline):
    """Infer sex from X and Y depth so the sex chromosomes are scored sanely.

    Without this a normal male is reported as X and Y at zero copies, because
    both sit near half the autosomal level and the copy-number model expects a
    diploid reference everywhere. Returns (label, expected_copies_for_X_and_Y).
    """
    def level(chrom):
        selector = bins["chrom"] == chrom
        if not numpy.any(selector):
            return float("nan")
        values = bins["adjusted"][selector]
        values = values[numpy.isfinite(values)]
        if values.shape[0] < 10:
            return float("nan")
        return float(numpy.median(values)) / baseline

    x_level, y_level = level("chrX"), level("chrY")
    if math.isfinite(y_level) and y_level > 0.25:
        return "XY", 1
    if math.isfinite(x_level) and x_level < 0.75:
        return "XY", 1
    return "XX", 2


def arm_of(chrom, start, end):
    """Arm label for a bin, or None if it straddles the centromere."""
    if chrom not in T2T_CENTROMERES:
        return chrom[3:]
    cen_start, cen_end = T2T_CENTROMERES[chrom]
    if end <= cen_start:
        return chrom[3:] + "p"
    if start >= cen_end:
        return chrom[3:] + "q"
    return None


def _arm_sort_key(label):
    """Order arms as 1p, 1q, 2p, ... 22, X, Y."""
    body = label.rstrip("pq")
    suffix = label[len(body):]
    if body == "X":
        index = 23
    elif body == "Y":
        index = 24
    else:
        try:
            index = int(body)
        except ValueError:
            index = 99
    return (index, {"p": 0, "q": 1, "": 2}.get(suffix, 3))


def iscn_string(arm_copy_numbers, sex="XX"):
    """Compose an ISCN-style seq karyotype string from arm-level copy numbers.

    p and q are merged when they carry the same copy number, matching the
    convention used by lrdk.py so that the two can be compared directly.

    Note on a defect inherited from that implementation and fixed here: the
    original collapses the diploid group to the literal "1-22" only when every
    autosome is diploid, but still enters that branch whenever X or Y is
    diploid. In any sample that is partly aneuploid with a diploid X — which is
    most female myeloma genomes — every diploid autosome was silently dropped
    from the output, leaving "(X)x2" alone. Groups are now always listed in
    full, with "1-22" used only as a genuine shorthand.
    """
    by_copy_number = defaultdict(list)
    for arm in arm_copy_numbers:
        if arm == "Y" and sex == "XX":
            continue  # a female sample has no Y to report
        by_copy_number[arm_copy_numbers[arm]].append(arm)

    # Merge p and q wherever both are present at the same copy number. Done by
    # inspection rather than against a fixed list of metacentrics, so it stays
    # correct if the centromere table changes.
    for copy_number in by_copy_number:
        arms = by_copy_number[copy_number]
        for chrom in [str(c) for c in range(1, 23)] + ["X", "Y"]:
            p, q = chrom + "p", chrom + "q"
            if p in arms and q in arms:
                arms.remove(p)
                arms.remove(q)
                arms.append(chrom)
        arms.sort(key=_arm_sort_key)

    autosome_labels = set(str(c) for c in range(1, 23))
    parts = []
    for copy_number in sorted(by_copy_number):
        arms = by_copy_number[copy_number]
        if not arms:
            continue
        if copy_number == 2 and autosome_labels.issubset(set(arms)):
            remainder = [a for a in arms if a not in autosome_labels]
            arms = ["1-22"] + remainder
        parts.append("(%s)x%d" % (", ".join(arms), copy_number))
    return "seq" + " ".join(parts)


def load_cytobands(path, with_stain=False):
    """Read a cytoband BED (chrom, start, end, name[, stain])."""
    bands = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom = canonical_chrom(fields[0])
            if chrom is None:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if with_stain:
                bands.append((chrom, start, end, fields[3],
                              fields[4] if len(fields) > 4 else "gneg"))
            else:
                bands.append((chrom, start, end, fields[3]))
    return bands


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_genome(bins, segments, output_path, sample_label, fit):
    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot
    except ImportError:
        sys.stderr.write("WARNING: matplotlib unavailable, skipping plot\n")
        return False

    figure, axis = pyplot.subplots(figsize=(20, 5), dpi=110)
    offset = 0
    tick_positions, tick_labels, separators = [], [], []
    offsets = {}

    for chrom in bins["order"]:
        indices = bins["chrom"] == chrom
        count = int(numpy.count_nonzero(indices))
        if count == 0:
            continue
        offsets[chrom] = offset
        x = offset + numpy.arange(count)
        axis.scatter(x, bins["log2"][indices], s=1.5, color="#4C72B0",
                     alpha=0.45, linewidths=0)
        tick_positions.append(offset + count / 2.0)
        tick_labels.append(chrom[3:])
        offset += count
        separators.append(offset)

    for segment in segments:
        chrom = segment["chrom"]
        if chrom not in offsets:
            continue
        x0 = offsets[chrom] + segment["bin_lo"]
        x1 = offsets[chrom] + segment["bin_hi"]
        colour = "#C44E52" if abs(segment["log2"]) >= 0.15 else "#55A868"
        axis.plot([x0, x1], [segment["log2"], segment["log2"]],
                  color=colour, linewidth=2.5, solid_capstyle="butt")

    for position in separators[:-1]:
        axis.axvline(position, color="0.8", linewidth=0.6)
    axis.axhline(0.0, color="0.3", linewidth=0.8, linestyle="--")

    axis.set_xticks(tick_positions)
    axis.set_xticklabels(tick_labels, fontsize=9)
    axis.set_xlim(0, offset)
    axis.set_ylim(-2.0, 2.0)
    axis.set_ylabel("log2 read-count ratio")
    axis.set_xlabel("Chromosome")
    axis.set_title("%s   purity %.2f   ploidy %.2f"
                   % (sample_label, fit["purity"], fit["ploidy"]))
    figure.tight_layout()
    figure.savefig(output_path)
    pyplot.close(figure)
    return True


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def load_fai(path):
    lengths = OrderedDict()
    contig_names = {}
    with open(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            chrom = canonical_chrom(fields[0])
            if chrom is None:
                continue
            lengths[chrom] = int(fields[1])
            contig_names[chrom] = fields[0]
    return lengths, contig_names


def main():
    parser = argparse.ArgumentParser(
        description="Digital karyotyping and CNA from nanopore adaptive "
                    "sampling BAM files (T2T-CHM13v2.0)")
    parser.add_argument("bam", help="Coordinate-sorted BAM aligned to T2T-CHM13v2.0")
    parser.add_argument("fai", help="FASTA index (.fai) for the alignment reference")
    parser.add_argument("prefix", help="Output path prefix")
    parser.add_argument("--exclude", action="append", default=[], metavar="BED",
                        help="BED of regions to exclude. Repeatable. Pass the "
                             "adaptive sampling panel BED here, plus repeat "
                             "and satellite annotation.")
    parser.add_argument("--reference", metavar="FASTA",
                        help="Reference FASTA for GC correction (optional but "
                             "recommended)")
    parser.add_argument("--cytobands", metavar="BED",
                        help="Cytoband BED for per-band reporting (optional)")
    parser.add_argument("--mask-stains", default="acen,stalk,gvar",
                        help="Cytoband stains to exclude when --cytobands is "
                             "given. These are the acrocentric short arms, "
                             "centromeres and heterochromatic variants, where "
                             "rDNA and satellite make read counts meaningless "
                             "(default: acen,stalk,gvar). Empty string to keep "
                             "them.")
    parser.add_argument("--min-band-bins", type=int, default=5,
                        help="Minimum bins for a cytoband row (default: 5)")
    parser.add_argument("--bin-size", type=int, default=200000,
                        help="Bin size in bp (default: 200000)")
    parser.add_argument("--mask-resolution", type=int, default=1000,
                        help="Mask granularity in bp (default: 1000)")
    parser.add_argument("--min-unmasked", type=float, default=0.25,
                        help="Minimum unmasked fraction to keep a bin "
                             "(default: 0.25)")
    parser.add_argument("--min-mapq", type=int, default=1,
                        help="Minimum mapping quality (default: 1)")
    parser.add_argument("--min-seg-bins", type=int, default=5,
                        help="Minimum bins per segment (default: 5)")
    parser.add_argument("--min-delta", type=float, default=0.07,
                        help="Minimum log2 difference between adjacent "
                             "segments. 0.07 corresponds to a single-copy "
                             "change in about 10 percent of cells; 0.10 to "
                             "about 14 percent (default: 0.07)")
    parser.add_argument("--seg-threshold", type=float, default=4.0,
                        help="t statistic threshold for a split (default: 4.0)")
    parser.add_argument("--baseline", choices=["mode", "median"], default="mode",
                        help="Reference level for the log2 ratio. 'mode' is "
                             "the dominant copy state and is the correct "
                             "choice for aneuploid genomes; 'median' is biased "
                             "upward when much of the genome is gained "
                             "(default: mode)")
    parser.add_argument("--purity", type=float, default=None,
                        help="Externally measured tumour purity, for example "
                             "post-sort flow plasma cell fraction as 0-1. "
                             "Fixes purity instead of fitting it.")
    parser.add_argument("--max-cn", type=int, default=8,
                        help="Highest integer copy state considered (default: 8)")
    parser.add_argument("--ploidy-tolerance", type=float, default=0.01,
                        help="Score slack within which the lowest ploidy is "
                             "preferred (default: 0.01)")
    parser.add_argument("--no-ploidy-parsimony", action="store_true",
                        help="Report the unconstrained best fit even when a "
                             "lower-ploidy solution fits comparably well")
    parser.add_argument("--no-centromere-mask", action="store_true",
                        help="Do not auto-exclude centromere spans")
    parser.add_argument("--threads", type=int, default=4,
                        help="BAM decompression threads (default: 4)")
    parser.add_argument("--read-limit", type=int, default=None,
                        help="Stop after this many BAM records (for testing)")
    parser.add_argument("--allow-truncated", action="store_true",
                        help="Tolerate a truncated BAM so the script can run "
                             "against an in-progress sequencing run. Disables "
                             "multithreaded decompression.")
    parser.add_argument("--sample", default=None,
                        help="Sample label for outputs (default: BAM basename)")
    parser.add_argument("--sex", choices=["auto", "XX", "XY"], default="auto",
                        help="Sample sex, for the expected copy number on X and "
                             "Y. Default 'auto' infers it from X and Y depth. "
                             "Supply it when known: a wrong inference turns a "
                             "normal male X into a whole-chromosome deletion. "
                             "A supplied value that disagrees with the depth "
                             "inference is reported as a warning, which is a "
                             "sample-identity check worth attending to.")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.bin_size % args.mask_resolution != 0:
        raise SystemExit("ERROR: --bin-size must be a multiple of "
                         "--mask-resolution")

    sample = args.sample or os.path.basename(args.bam).split(".")[0]
    output_dir = os.path.dirname(os.path.abspath(args.prefix))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    # --- reference and bins -------------------------------------------------
    sys.stderr.write("Reading reference index...\n")
    lengths, fasta_contigs = load_fai(args.fai)
    if not lengths:
        raise SystemExit("ERROR: no recognisable contigs in %s" % args.fai)
    genome = Genome(lengths, args.bin_size, args.mask_resolution)
    genome.fasta_contig = fasta_contigs

    # --- mask ---------------------------------------------------------------
    sys.stderr.write("Building exclusion mask...\n")
    if not args.no_centromere_mask:
        genome.mask_centromeres()
    for bed in args.exclude:
        used = genome.add_mask_bed(bed)
        sys.stderr.write("  %s: %d intervals\n" % (bed, used))
    if not args.exclude:
        sys.stderr.write(
            "WARNING: no --exclude BED supplied. The adaptive sampling panel "
            "will dominate the signal and the result will be meaningless.\n")

    # Acrocentric short arms, centromeres and heterochromatic variants carry
    # rDNA arrays and satellite. In T2T they are assembled but not uniquely
    # mappable, so read counts there are meaningless: on the first banded run
    # 13p, 14p, 15p, 21p, 22p and chrY came back at log2 -2 to -8, which is
    # "no reads" rather than "deleted". Masking them by Giemsa stain is
    # precise and needs no extra annotation beyond the cytoband file.
    if args.cytobands and args.mask_stains:
        wanted = {t.strip().lower() for t in args.mask_stains.split(",")
                  if t.strip()}
        masked_bands = 0
        for chrom, start, end, _name, stain in load_cytobands(args.cytobands,
                                                              with_stain=True):
            if stain.lower() not in wanted or chrom not in genome.mask:
                continue
            lo = start // genome.mask_resolution
            hi = min(end // genome.mask_resolution + 1,
                     genome.mask[chrom].shape[0])
            genome.mask[chrom][lo:hi] = True
            masked_bands += 1
        sys.stderr.write("  masked %d cytobands stained %s\n"
                         % (masked_bands, ", ".join(sorted(wanted))))

    fractions = genome.unmasked_fraction()
    total_excluded = 1.0 - numpy.mean(numpy.concatenate(
        [fractions[c] for c in fractions]))
    sys.stderr.write("  genome excluded: %.1f%%\n" % (100.0 * total_excluded))

    # --- counting -----------------------------------------------------------
    sys.stderr.write("Counting reads...\n")
    stats = count_reads(args.bam, genome, args.min_mapq, args.threads,
                        args.read_limit, args.verbose, args.allow_truncated)
    sys.stderr.write("  %d records, %d counted, %d in masked regions\n"
                     % (stats["total"], stats["counted"], stats["masked"]))
    if stats["counted"] == 0:
        raise SystemExit("ERROR: no reads counted. Check the reference, the "
                         "mask and --min-mapq.")

    # --- per-bin normalisation ---------------------------------------------
    sys.stderr.write("Normalising bins...\n")
    fasta = None
    if args.reference:
        fasta = pysam.FastaFile(args.reference)
        missing = [genome.fasta_contig[c] for c in genome.lengths
                   if genome.fasta_contig[c] not in set(fasta.references)]
        if missing:
            fasta.close()
            raise SystemExit(
                "ERROR: --reference does not contain contigs named in the .fai "
                "(for example %s). Pass the FASTA that this .fai indexes."
                % missing[0])

    chrom_column, start_column, end_column = [], [], []
    raw_column, frac_column, gc_column, adjusted_column = [], [], [], []

    for chrom in genome.lengths:
        n_bins = genome.n_bins[chrom]
        counts = genome.counts[chrom].astype(numpy.float64)
        fraction = fractions[chrom]
        keep = fraction >= args.min_unmasked
        if not numpy.any(keep):
            continue
        # Scale each bin up to what it would have been at full width.
        adjusted = numpy.where(keep, counts / numpy.maximum(fraction, 1e-9), 0.0)
        if fasta is not None:
            gc = gc_content(fasta, genome, chrom, n_bins)
        else:
            gc = numpy.full(n_bins, numpy.nan)

        indices = numpy.nonzero(keep)[0]
        chrom_column.extend([chrom] * indices.shape[0])
        start_column.append(indices * genome.bin_size)
        end_column.append(numpy.minimum((indices + 1) * genome.bin_size,
                                        genome.lengths[chrom]))
        raw_column.append(counts[indices])
        frac_column.append(fraction[indices])
        gc_column.append(gc[indices])
        adjusted_column.append(adjusted[indices])

    if fasta is not None:
        fasta.close()

    bins = {
        "chrom": numpy.array(chrom_column),
        "start": numpy.concatenate(start_column).astype(numpy.int64),
        "end": numpy.concatenate(end_column).astype(numpy.int64),
        "raw": numpy.concatenate(raw_column),
        "unmasked": numpy.concatenate(frac_column),
        "gc": numpy.concatenate(gc_column),
        "adjusted": numpy.concatenate(adjusted_column),
    }
    bins["order"] = [c for c in genome.lengths
                     if numpy.any(bins["chrom"] == c)]

    # GC correction across the whole genome at once.
    gc_applied = False
    if fasta is not None:
        bins["adjusted"], gc_applied = correct_gc(bins["adjusted"], bins["gc"])
    sys.stderr.write("  GC correction: %s\n"
                     % ("applied" if gc_applied else "not applied"))

    autosomal = numpy.array([c not in ("chrX", "chrY") for c in bins["chrom"]])
    positive = bins["adjusted"] > 0
    reference_values = bins["adjusted"][autosomal & positive]
    if args.baseline == "median":
        baseline = float(numpy.median(reference_values))
        baseline_method = "median"
        median_baseline = baseline
    else:
        # Estimate per-bin noise first so the mode finder can match its
        # bandwidth to it. Successive differencing is scale invariant, so this
        # works on the un-baselined values.
        sigma_hint = robust_sd_from_diffs(numpy.log2(reference_values))
        baseline, baseline_method, median_baseline = modal_baseline(
            reference_values, sigma_hint)
    if not (baseline > 0):
        raise SystemExit("ERROR: autosomal baseline is zero.")
    baseline_shift = (math.log2(median_baseline / baseline)
                      if baseline > 0 and median_baseline > 0 else float("nan"))
    sys.stderr.write("  baseline %.1f reads/bin (%s); median estimate %.1f "
                     "(%+.3f log2)\n"
                     % (baseline, baseline_method, median_baseline, baseline_shift))
    if math.isfinite(baseline_shift) and abs(baseline_shift) > 0.03:
        sys.stderr.write("  note: modal and median baselines differ, which "
                         "indicates a substantially aneuploid genome\n")
    with numpy.errstate(divide="ignore", invalid="ignore"):
        bins["log2"] = numpy.log2(bins["adjusted"] / baseline)
    bins["log2"][~numpy.isfinite(bins["log2"])] = numpy.nan

    sex_label, sex_reference_copies = infer_sex(bins, baseline)
    sex_source = "inferred"
    if args.sex != "auto":
        # Keep the inference and compare. A disagreement means the sequenced
        # material is not the sex the record says it is, which is a sample
        # identity problem rather than a copy-number one.
        if args.sex != sex_label:
            sys.stderr.write(
                "WARNING: supplied sex %s disagrees with the sex inferred from "
                "depth (%s) for sample %s. The supplied value is being used. "
                "Check that this tube holds the material it is labelled with "
                "before interpreting any copy number from it.\n"
                % (args.sex, sex_label, args.sample or args.prefix))
        sex_label = args.sex
        sex_reference_copies = 1 if args.sex == "XY" else 2
        sex_source = "supplied"
    sys.stderr.write("  inferred sex: %s\n" % sex_label)

    median_reads = float(numpy.median(bins["raw"][autosomal & positive]))
    sigma = robust_sd_from_diffs(bins["log2"][autosomal])
    poisson_sigma = (1.4427 / math.sqrt(median_reads)) if median_reads > 0 else float("nan")
    overdispersion = sigma / poisson_sigma if poisson_sigma > 0 else float("nan")
    sys.stderr.write("  median reads/bin %.0f, noise sigma %.4f "
                     "(Poisson floor %.4f, overdispersion %.2fx)\n"
                     % (median_reads, sigma, poisson_sigma, overdispersion))

    # --- segmentation -------------------------------------------------------
    sys.stderr.write("Segmenting...\n")
    segments = []
    for chrom in bins["order"]:
        selector = bins["chrom"] == chrom
        chrom_indices = numpy.nonzero(selector)[0]
        # Map a global bin index to its position within this chromosome, used
        # only for plotting coordinates.
        local_index = {int(g): k for k, g in enumerate(chrom_indices)}
        arms = numpy.array([arm_of(chrom, bins["start"][i], bins["end"][i])
                            for i in chrom_indices])
        for arm in [a for a in OrderedDict.fromkeys(arms) if a is not None]:
            arm_positions = chrom_indices[arms == arm]
            values = bins["log2"][arm_positions]
            finite = numpy.isfinite(values)
            if numpy.count_nonzero(finite) < args.min_seg_bins:
                continue
            filler = float(numpy.nanmedian(values))
            if not math.isfinite(filler):
                continue
            values = numpy.where(finite, values, filler)
            boundaries = segment_values(values, sigma, args.min_seg_bins,
                                        args.seg_threshold)
            boundaries = merge_segments(values, boundaries, args.min_delta)
            for lo, hi in boundaries:
                positions = arm_positions[lo:hi]
                window = bins["log2"][positions]
                if not numpy.any(numpy.isfinite(window)):
                    continue
                median_log2 = float(numpy.nanmedian(window))
                n_bins = int(positions.shape[0])
                standard_error = sigma / math.sqrt(n_bins)
                segments.append({
                    "chrom": chrom,
                    "arm": arm,
                    "start": int(bins["start"][positions[0]]),
                    "end": int(bins["end"][positions[-1]]),
                    "bin_lo": local_index[int(positions[0])],
                    "bin_hi": local_index[int(positions[-1])] + 1,
                    "n_bins": n_bins,
                    "log2": median_log2,
                    "ci_low": median_log2 - 1.96 * standard_error,
                    "ci_high": median_log2 + 1.96 * standard_error,
                    "reads": float(numpy.sum(bins["raw"][positions])),
                    "indices": positions,
                })

    sys.stderr.write("  %d segments\n" % len(segments))

    # --- purity and ploidy --------------------------------------------------
    sys.stderr.write("Fitting absolute copy number...\n")
    # The fit assumes a diploid germline reference, so the sex chromosomes are
    # held out of it; they are still assigned copy numbers afterwards.
    finite_segments = [s for s in segments
                       if math.isfinite(s["log2"])
                       and s["chrom"] not in ("chrX", "chrY")]
    if not finite_segments:
        raise SystemExit("ERROR: no segments with finite log2 ratios.")
    seg_log2 = numpy.array([s["log2"] for s in finite_segments])
    seg_weight = numpy.array([float(s["n_bins"]) for s in finite_segments])
    fit = fit_purity_ploidy(seg_log2, seg_weight, args.purity, args.max_cn,
                            tolerance=args.ploidy_tolerance,
                            parsimony=not args.no_ploidy_parsimony)
    sys.stderr.write("  purity %.3f, ploidy %.2f, score %.4f, "
                     "ploidy margin %.4f%s\n"
                     % (fit["purity"], fit["ploidy"], fit["score"],
                        fit["margin"],
                        " (purity fixed by user)" if fit["fixed_purity"] else ""))
    sys.stderr.write("  genome not explained by a single clone: %.1f%%\n"
                     % (100.0 * fit["unexplained_fraction"]))
    # The score is a length-weighted distance in units of the copy-state gap,
    # so it is bounded at 0.5 and anything above about 0.08 means the segments
    # are not landing on integer states.
    if fit["score"] > 0.08:
        sys.stderr.write(
            "WARNING: poor fit (score %.4f). Segments are not resolving to "
            "integer copy states. %s\n"
            % (fit["score"],
               "Check that --purity is the post-sort flow purity of the "
               "sorted cells, not the marrow plasma cell percentage."
               if fit["fixed_purity"] else
               "Consider supplying --purity from flow cytometry."))
    if fit["ambiguous"]:
        sys.stderr.write(
            "WARNING: purity and ploidy are not identifiable from depth alone "
            "on this sample. Reported (%.2f, %.2f) was chosen by ploidy "
            "parsimony; the unconstrained best fit was (%.2f, %.2f). Supply "
            "--purity from post-sort flow cytometry to resolve this.\n"
            % (fit["purity"], fit["ploidy"], fit["global_purity"],
               fit["global_ploidy"]))

    for segment in segments:
        copy_number, residual = assign_copy_number(
            segment["log2"], fit["purity"], fit["ploidy"], args.max_cn,
            sex_reference_copies if segment["chrom"] in ("chrX", "chrY") else 2)
        # chrY is not scored below whole-chromosome level. Yq11.222, q11.223
        # and q11.23 are ampliconic and palindromic, and read 1 to 3 log2
        # below the true single-copy level because reads do not map there
        # uniquely. Loss of Y is still called on the arm row.
        kind, fraction = classify_event(
            copy_number, segment["log2"],
            None if segment["chrom"] == "chrY"
            else event_expected_copies(segment["chrom"], sex_label))
        segment["cn"] = copy_number
        segment["cn_residual"] = residual
        segment["event"] = kind
        segment["cell_fraction"] = fraction
        segment["detect_limit"] = detection_limit(sigma, segment["n_bins"])

    # --- arm summary --------------------------------------------------------
    arm_rows = []
    arm_copy_numbers = OrderedDict()
    for chrom in bins["order"]:
        selector = bins["chrom"] == chrom
        chrom_indices = numpy.nonzero(selector)[0]
        arms = numpy.array([arm_of(chrom, bins["start"][i], bins["end"][i])
                            for i in chrom_indices])
        for arm in [a for a in OrderedDict.fromkeys(arms) if a is not None]:
            positions = chrom_indices[arms == arm]
            values = bins["log2"][positions]
            values = values[numpy.isfinite(values)]
            if values.shape[0] == 0:
                continue
            median_log2 = float(numpy.median(values))
            copy_number, _residual = assign_copy_number(
                median_log2, fit["purity"], fit["ploidy"], args.max_cn,
                sex_reference_copies if chrom in ("chrX", "chrY") else 2)
            # No special case for the sex chromosomes: event_expected_copies
            # supplies the expectation, so a normal male X at cn 1 returns
            # "none" for the same reason a normal autosome at cn 2 does, and
            # chrY in a female is not scored at all.
            kind, fraction = classify_event(
                copy_number, median_log2,
                event_expected_copies(chrom, sex_label))
            arm_copy_numbers[arm] = copy_number
            arm_rows.append({
                "arm": arm, "chrom": chrom, "n_bins": int(values.shape[0]),
                "log2": median_log2, "cn": copy_number, "event": kind,
                "cell_fraction": fraction,
                "detect_limit": detection_limit(sigma, int(values.shape[0])),
            })

    # --- hyperdiploidy ------------------------------------------------------
    # Gains are judged against the copy number the modal baseline maps to, not
    # against a hardcoded 3. On the first real sample a ploidy-3 fit made every
    # baseline chromosome satisfy cn >= 3 and the trisomy list filled up with
    # chromosomes sitting at the baseline.
    baseline_cn, _residual = assign_copy_number(
        0.0, fit["purity"], fit["ploidy"], args.max_cn)
    trisomic = []
    for row in arm_rows:
        if row["chrom"] in HYPERDIPLOID_CHROMS and row["cn"] > baseline_cn:
            label = row["chrom"][3:]
            if label not in trisomic:
                trisomic.append(label)
    hyperdiploid = len(trisomic) >= 2

    # --- writing ------------------------------------------------------------
    sys.stderr.write("Writing outputs...\n")

    with open(args.prefix + ".bins.tsv", "w") as handle:
        handle.write("chromosome\tstart\tend\tn_reads\tunmasked_fraction\t"
                     "gc\tadjusted\tlog2_ratio\n")
        for i in range(bins["chrom"].shape[0]):
            handle.write("%s\t%d\t%d\t%d\t%.4f\t%s\t%.2f\t%s\n" % (
                bins["chrom"][i], bins["start"][i], bins["end"][i],
                int(bins["raw"][i]), bins["unmasked"][i],
                ("%.4f" % bins["gc"][i]) if numpy.isfinite(bins["gc"][i]) else "NA",
                bins["adjusted"][i],
                ("%.4f" % bins["log2"][i]) if numpy.isfinite(bins["log2"][i]) else "NA"))

    with open(args.prefix + ".segments.tsv", "w") as handle:
        handle.write("chromosome\tstart\tend\tarm\tn_bins\tsize_mb\treads\t"
                     "log2_ratio\tci_low\tci_high\tcopy_number\tcn_residual\t"
                     "event\tcell_fraction\tdetection_limit\n")
        for segment in segments:
            handle.write("%s\t%d\t%d\t%s\t%d\t%.2f\t%d\t%.4f\t%.4f\t%.4f\t"
                         "%d\t%.4f\t%s\t%s\t%s\n" % (
                segment["chrom"], segment["start"], segment["end"],
                segment["arm"], segment["n_bins"],
                (segment["end"] - segment["start"]) / 1e6, int(segment["reads"]),
                segment["log2"], segment["ci_low"], segment["ci_high"],
                segment["cn"], segment["cn_residual"], segment["event"],
                ("%.3f" % segment["cell_fraction"])
                if math.isfinite(segment["cell_fraction"]) else "NA",
                ("%.3f" % segment["detect_limit"])
                if math.isfinite(segment["detect_limit"]) else "NA"))

    with open(args.prefix + ".arms.tsv", "w") as handle:
        handle.write("arm\tchromosome\tn_bins\tlog2_ratio\tcopy_number\t"
                     "event\tcell_fraction\tdetection_limit\n")
        for row in arm_rows:
            handle.write("%s\t%s\t%d\t%.4f\t%d\t%s\t%s\t%s\n" % (
                row["arm"], row["chrom"], row["n_bins"], row["log2"], row["cn"],
                row["event"],
                ("%.3f" % row["cell_fraction"])
                if math.isfinite(row["cell_fraction"]) else "NA",
                ("%.3f" % row["detect_limit"])
                if math.isfinite(row["detect_limit"]) else "NA"))

    cytoband_rows = []
    if args.cytobands:
        bands = load_cytobands(args.cytobands)
        with open(args.prefix + ".cytobands.tsv", "w") as handle:
            handle.write("cytoband\tchromosome\tstart\tend\tn_bins\t"
                         "log2_ratio\tcopy_number\timplied_copies\tevent\t"
                         "cell_fraction\tdetection_limit\n")
            for chrom, start, end, name in bands:
                selector = ((bins["chrom"] == chrom) &
                            (bins["start"] < end) & (bins["end"] > start))
                if chrom == "chrY" and sex_label == "XX":
                    continue
                values = bins["log2"][selector]
                values = values[numpy.isfinite(values)]
                # A band resting on one or two bins has a detection limit worse
                # than 0.2, so it cannot support a call.
                if values.shape[0] < args.min_band_bins:
                    continue
                median_log2 = float(numpy.median(values))
                copy_number, _residual = assign_copy_number(
                    median_log2, fit["purity"], fit["ploidy"], args.max_cn,
                    sex_reference_copies if chrom in ("chrX", "chrY") else 2)
                # See the segment loop: chrY is scored at arm level only.
                kind, fraction = classify_event(
                    copy_number, median_log2,
                    None if chrom == "chrY"
                    else event_expected_copies(chrom, sex_label))
                label = chrom[3:] + name
                row = {"cytoband": label, "log2": median_log2,
                       "cn": copy_number, "event": kind,
                       "cell_fraction": fraction,
                       "detect_limit": detection_limit(sigma, values.shape[0])}
                cytoband_rows.append(row)
                # Decimal copies is always defined; cell_fraction is only
                # meaningful for a single-copy change on a diploid background
                # and is NA above log2 0.585, where more than one extra copy is
                # required. Reporting both avoids a column of bare NA on every
                # trisomy that has drifted toward tetrasomy.
                reference = 1 if chrom in ("chrX", "chrY") and sex_label == "XY" else 2
                implied = ((2.0 ** median_log2)
                           * (fit["purity"] * fit["ploidy"]
                              + (1 - fit["purity"]) * 2.0)
                           - reference * (1 - fit["purity"])) / fit["purity"]
                handle.write("%s\t%s\t%d\t%d\t%d\t%.4f\t%d\t%.2f\t%s\t%s\t%s\n" % (
                    label, chrom, start, end, values.shape[0], median_log2,
                    copy_number, implied, kind,
                    ("%.3f" % fraction) if math.isfinite(fraction) else "NA",
                    ("%.3f" % row["detect_limit"])
                    if math.isfinite(row["detect_limit"]) else "NA"))

    summary = {
        "sample": sample,
        "bam": os.path.abspath(args.bam),
        "parameters": {
            "bin_size": args.bin_size,
            "mask_resolution": args.mask_resolution,
            "min_unmasked": args.min_unmasked,
            "min_mapq": args.min_mapq,
            "min_seg_bins": args.min_seg_bins,
            "min_delta": args.min_delta,
            "seg_threshold": args.seg_threshold,
            "max_cn": args.max_cn,
            "exclude_beds": [os.path.abspath(b) for b in args.exclude],
            "gc_correction": bool(gc_applied),
            "centromere_mask": not args.no_centromere_mask,
        },
        "qc": {
            "records_read": stats["total"],
            "reads_counted": stats["counted"],
            "reads_in_masked_regions": stats["masked"],
            "reads_below_mapq": stats["low_mapq"],
            "secondary_or_supplementary": stats["secondary_supp"],
            "bins_retained": int(bins["chrom"].shape[0]),
            "genome_excluded_fraction": round(float(total_excluded), 4),
            "median_reads_per_bin": round(median_reads, 1),
            "baseline_reads_per_bin": round(float(baseline), 1),
            "baseline_method": baseline_method,
            "baseline_vs_median_log2": (round(baseline_shift, 4)
                                        if math.isfinite(baseline_shift) else None),
            "noise_sigma_log2": round(sigma, 5),
            "poisson_floor_log2": round(poisson_sigma, 5),
            "overdispersion": round(overdispersion, 3),
            "detection_limit_5mb": round(
                detection_limit(sigma, max(1, 5000000 // args.bin_size)), 4),
            "detection_limit_arm": round(detection_limit(sigma, 200), 4),
        },
        "fit": {
            "purity": round(fit["purity"], 4),
            "ploidy": round(fit["ploidy"], 3),
            "score": round(fit["score"], 5),
            "ploidy_margin": (round(fit["margin"], 5)
                              if math.isfinite(fit["margin"]) else None),
            "unexplained_fraction": round(fit["unexplained_fraction"], 4),
            "ambiguous": bool(fit["ambiguous"]),
            "unconstrained_best": {"purity": round(fit["global_purity"], 4),
                                   "ploidy": round(fit["global_ploidy"], 3),
                                   "score": round(fit["global_score"], 5)},
            "purity_supplied": fit["fixed_purity"],
            "poor_fit": bool(fit["score"] > 0.08),
        },
        "inferred_sex": sex_label,
        "hyperdiploidy": {
            "called": hyperdiploid,
            "baseline_copy_number": baseline_cn,
            "trisomic_chromosomes": trisomic,
            "panel": [c[3:] for c in HYPERDIPLOID_CHROMS],
        },
        "arms": {row["arm"]: {
            "cn": row["cn"],
            "log2": round(row["log2"], 4),
            "event": row["event"],
            "cell_fraction": (round(row["cell_fraction"], 3)
                              if math.isfinite(row["cell_fraction"]) else None),
        } for row in arm_rows},
        "ISCN_karyotype": iscn_string(arm_copy_numbers, sex_label),
        "n_segments": len(segments),
        "altered_segments": [
            {"region": "%s:%d-%d" % (s["chrom"], s["start"], s["end"]),
             "size_mb": round((s["end"] - s["start"]) / 1e6, 2),
             "log2": round(s["log2"], 4),
             "cn": s["cn"],
             "event": s["event"],
             "cell_fraction": (round(s["cell_fraction"], 3)
                               if math.isfinite(s["cell_fraction"]) else None)}
            # Filtered on the event, not on abs(log2) >= min_delta. A region
            # sitting at its fitted copy number is not altered whatever its
            # log2 is, and the old threshold put two thirds of a near-diploid
            # genome into a list named "altered".
            for s in segments if s["event"] != "none"
        ],
    }
    if cytoband_rows:
        summary["cytobands_altered"] = [
            {"band": r["cytoband"], "cn": r["cn"], "log2": round(r["log2"], 4),
             "event": r["event"],
             "cell_fraction": (round(r["cell_fraction"], 3)
                               if math.isfinite(r["cell_fraction"]) else None)}
            for r in cytoband_rows if r["event"] != "none"
        ]

    with open(args.prefix + ".karyotype.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    plot_genome(bins, segments, args.prefix + ".karyotype.png", sample, fit)

    sys.stderr.write("\n%s\n" % summary["ISCN_karyotype"])
    sys.stderr.write("Hyperdiploid: %s%s\n" % (
        hyperdiploid,
        (" (trisomies %s)" % ", ".join(trisomic)) if trisomic else ""))
    sys.stderr.write("Detection limit, 5 Mb event: %.1f%% of cells\n"
                     % (100.0 * summary["qc"]["detection_limit_5mb"]))
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
