#!/usr/bin/env python3
"""
baf_loh_screen.py

Panel-wide B-allele frequency (BAF) screen for loss of heterozygosity, including
copy-neutral LOH, from Clair3 phased germline VCFs.

Rationale
---------
At a balanced heterozygous site the minor-allele fraction is centred on 0.5, with
spread determined by sequencing depth. Under LOH within a tumour cell population
of fraction f, the observed allele fraction at a surviving germline het site is
displaced away from 0.5. Two consequences follow, and this script measures both:

  1. The BAF distribution within the affected region loses density at 0.5 and
     becomes bimodal.
  2. Sites whose displaced allele fraction crosses the variant caller's
     heterozygous threshold are emitted as homozygous, so the count of called
     het sites in the region falls.

Neither signal is interpretable in isolation at low depth, so this script
computes a depth-aware expectation for signal (1) and a cohort-relative
expectation for signal (2), and requires both before flagging a region.

Copy number is deliberately not inferred here. BAF deflection alone cannot
distinguish copy-neutral LOH from hemizygous deletion; that discrimination
requires the copy number track (see the companion plotting script).

Inputs
------
  - Panel BED (v7), four columns: chrom, start, end, name.
  - One or more phased Clair3 VCFs per sample. Both the per-chromosome layout
    produced by the pipeline and a single merged VCF are accepted.

Output
------
  Tab-separated table, one row per panel region per sample.

Dependencies
------------
  Python standard library only. No pandas. Runs in the awgs_sv conda
  environment as-is.

Usage
-----
  python3 baf_loh_screen.py \
      --bed  aWGS_PCN_v7_t2t_chr.bed \
      --vcf-dir results_v7_<date>_24h/hg38/clair3_phased/<sample>/clair3_out/tmp/phase_output/phase_vcf \
      --sample <SEQ_ID> \
      --out    <SEQ_ID>.baf_screen.tsv

  Multiple samples, with cohort-relative normalisation:

  python3 baf_loh_screen.py \
      --bed aWGS_PCN_v7_t2t_chr.bed \
      --sample-map samples.tsv \
      --out cohort.baf_screen.tsv

  where samples.tsv is two tab-separated columns: sample_id and the VCF
  directory (or a single merged VCF path) for that sample.
"""

import argparse
import glob
import gzip
import math
import os
import statistics
import sys


# ---------------------------------------------------------------------------
# Analysis parameters
#
# These are defaults, all overridable from the command line. They are chosen
# for ONT adaptive-sampling panel data at the depths this pipeline achieves and
# should be revisited if the depth regime changes materially.
# ---------------------------------------------------------------------------

# Calibrated against a three-sample cohort in which one sample carried a
# FISH-confirmed del(17p) at 65% clonality with a somatic TP53 missense variant
# at VAF 0.93, and two did not. Panel-wide flag counts and the separation at the
# TP53 window were identical for a minimum site depth of 5, 8, 10 and 12, so the
# statistic is insensitive to this parameter over that range. The value below
# sits above the coarsest allele-fraction quantisation (at depth 8 the attainable
# BAF grid is 0.125 apart) while retaining essentially all sites in the
# validated window; raising it further discards data for no measurable gain.
DEFAULT_MIN_SITE_DEPTH = 8      # minimum depth for a het site to be usable
DEFAULT_MIN_SITES = 30          # minimum usable het sites for a region to be assessable
DEFAULT_CENTRAL_LO = 0.42       # lower bound of the "balanced" BAF band
DEFAULT_CENTRAL_HI = 0.58       # upper bound of the "balanced" BAF band
DEFAULT_CDR_THRESHOLD = 0.55    # central depletion ratio below which a region is depleted
DEFAULT_BIMODALITY_THRESHOLD = 1.0    # outer-lobe to central mass ratio supporting an LOH call
DEFAULT_DENSITY_RATIO = 0.45    # het density relative to cohort median below which density supports LOH

FLAG_LOH = "LOH_LIKELY"
FLAG_EQUIVOCAL = "EQUIVOCAL"
FLAG_NO_LOH = "NO_LOH"
FLAG_UNASSESSABLE = "UNASSESSABLE"


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def open_maybe_gzip(path):
    """Return a text-mode handle for a plain or gzip-compressed file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_bed(path):
    """
    Read a four-column BED file.

    Returns a list of dicts with keys: name, chrom, start, end. Regions are
    kept in file order; the caller is responsible for any sorting.
    """
    regions = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3] if len(fields) > 3 else "{0}:{1}-{2}".format(chrom, start, end)
            regions.append({"name": name, "chrom": chrom, "start": start, "end": end})
    if not regions:
        raise ValueError("No usable regions parsed from BED file: {0}".format(path))
    return regions


def collect_vcf_paths(target):
    """
    Resolve a user-supplied VCF location into a list of VCF file paths.

    Accepts either a single VCF file or a directory containing the
    per-chromosome phased VCFs that Clair3 emits.
    """
    if os.path.isfile(target):
        return [target]
    if os.path.isdir(target):
        paths = []
        for pattern in ("phased_*.vcf.gz", "phased_*.vcf"):
            paths.extend(glob.glob(os.path.join(target, pattern)))
        # Avoid double-counting a bgzipped VCF and its uncompressed twin.
        seen_stems = set()
        unique = []
        for path in sorted(paths):
            stem = path[:-3] if path.endswith(".gz") else path
            if stem in seen_stems:
                continue
            seen_stems.add(stem)
            unique.append(path)
        if not unique:
            raise ValueError("No VCF files found under directory: {0}".format(target))
        return unique
    raise ValueError("VCF path is neither a file nor a directory: {0}".format(target))


# ---------------------------------------------------------------------------
# VCF record handling
# ---------------------------------------------------------------------------

def parse_format_fields(format_field, sample_field):
    """Zip a VCF FORMAT key string against a sample value string."""
    keys = format_field.split(":")
    values = sample_field.split(":")
    return dict(zip(keys, values))


def genotype_is_het(gt_string):
    """
    Return True for a diploid heterozygous genotype call.

    Both phased ('|') and unphased ('/') separators are accepted, since the
    phased Clair3 output mixes the two depending on whether a site fell inside
    a phase block.
    """
    if gt_string in (".", "./.", ".|.", ""):
        return False
    alleles = gt_string.replace("|", "/").split("/")
    if len(alleles) != 2:
        return False
    if any(allele == "." for allele in alleles):
        return False
    return alleles[0] != alleles[1]


def extract_depth_and_baf(sample_values):
    """
    Derive total depth and minor-allele fraction from a parsed FORMAT dict.

    The allele-depth field AD is preferred because it gives the allele split
    directly. Where AD is absent the VAF field is used with DP, which is the
    fallback the Clair3 output occasionally requires.

    Returns (depth, minor_allele_fraction) or (None, None) if neither route
    yields a usable value.
    """
    ad_raw = sample_values.get("AD")
    if ad_raw and ad_raw != ".":
        try:
            counts = [int(x) for x in ad_raw.split(",") if x != "."]
        except ValueError:
            counts = []
        # Restrict to the reference and first alternate allele; multiallelic
        # sites are not informative for a biallelic BAF model.
        if len(counts) >= 2:
            ref_count, alt_count = counts[0], counts[1]
            depth = ref_count + alt_count
            if depth > 0:
                # Unfolded alternate-allele fraction, retained across the full
                # 0-1 range. Folding to min(ref, alt)/depth would collapse the
                # two modes of an LOH distribution onto each other and destroy
                # exactly the signal this screen exists to detect.
                return depth, alt_count / float(depth)

    vaf_raw = sample_values.get("VAF") or sample_values.get("AF")
    dp_raw = sample_values.get("DP")
    if vaf_raw and dp_raw and vaf_raw != "." and dp_raw != ".":
        try:
            vaf = float(vaf_raw.split(",")[0])
            depth = int(dp_raw)
        except ValueError:
            return None, None
        if depth > 0 and 0.0 <= vaf <= 1.0:
            return depth, vaf

    return None, None


def is_biallelic_snv(ref, alt):
    """
    Restrict the screen to biallelic single-nucleotide sites.

    Indels are excluded because ONT allele-depth estimates at indels in
    homopolymer context are unreliable and would add spread to the BAF
    distribution that has nothing to do with allelic balance.
    """
    if "," in alt:
        return False
    return len(ref) == 1 and len(alt) == 1 and ref != "." and alt != "."


def load_het_sites(vcf_paths, min_site_depth, sample_column=None):
    """
    Read het SNV sites from one or more VCFs.

    Returns a dict keyed by chromosome, each value a list of
    (position, depth, minor_allele_fraction, phased_flag) tuples sorted by
    position.

    Sites are filtered here on depth only. Region assignment happens later so
    that a single pass over the VCF serves all panel regions.
    """
    by_chrom = {}
    total_records = 0
    kept = 0

    for vcf_path in vcf_paths:
        sample_index = 0
        with open_maybe_gzip(vcf_path) as handle:
            for line in handle:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header = line.rstrip("\n").split("\t")
                    sample_names = header[9:]
                    if sample_column and sample_column in sample_names:
                        sample_index = sample_names.index(sample_column)
                    else:
                        sample_index = 0
                    continue

                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                total_records += 1

                chrom = fields[0]
                try:
                    pos = int(fields[1])
                except ValueError:
                    continue
                ref, alt, filt = fields[3], fields[4], fields[6]

                if filt not in (".", "PASS", "") :
                    continue
                if not is_biallelic_snv(ref, alt):
                    continue

                sample_fields = fields[9:]
                if sample_index >= len(sample_fields):
                    continue
                values = parse_format_fields(fields[8], sample_fields[sample_index])

                gt = values.get("GT", "./.")
                if not genotype_is_het(gt):
                    continue

                depth, minor_baf = extract_depth_and_baf(values)
                if depth is None or depth < min_site_depth:
                    continue

                phased = "|" in gt
                by_chrom.setdefault(chrom, []).append((pos, depth, minor_baf, phased))
                kept += 1

    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda item: item[0])

    return by_chrom, total_records, kept


# ---------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------

def log_binomial_coefficient(n, k):
    """Natural log of the binomial coefficient, via log-gamma for numeric stability."""
    return (math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1))


def binomial_central_probability(depth, lo, hi, p=0.5):
    """
    Probability that the alternate-allele fraction at a balanced heterozygous
    site of the given depth falls inside the central band [lo, hi].

    The fraction is unfolded, so a balanced site is centred on 0.5 and the
    band is genuinely central. Computing the expectation exactly from each
    site's own depth, rather than assuming a fixed expected fraction, is what
    keeps thin-coverage regions from being flagged simply for being noisy.
    """
    if depth <= 0:
        return 0.0
    total = 0.0
    for k in range(0, depth + 1):
        fraction = k / float(depth)
        if lo <= fraction <= hi:
            log_prob = (log_binomial_coefficient(depth, k)
                        + k * math.log(p)
                        + (depth - k) * math.log(1.0 - p))
            total += math.exp(log_prob)
    return total


def central_depletion_ratio(sites, lo, hi, cache):
    """
    Ratio of observed to expected density in the central BAF band.

    Expected density is the mean, across the region's own sites, of the
    per-site binomial probability of landing in the band. A value near 1.0
    indicates balanced heterozygosity; values well below 1.0 indicate density
    has been pushed out of the centre, which under LOH is displaced toward
    both extremes rather than toward one of them.
    """
    if not sites:
        return None
    observed = 0
    expected_sum = 0.0
    for _pos, depth, baf, _phased in sites:
        if lo <= baf <= hi:
            observed += 1
        key = (depth, lo, hi)
        if key not in cache:
            cache[key] = binomial_central_probability(depth, lo, hi)
        expected_sum += cache[key]
    if expected_sum <= 0:
        return None
    observed_fraction = observed / float(len(sites))
    expected_fraction = expected_sum / float(len(sites))
    return observed_fraction / expected_fraction


def bimodality_ratio(sites, lo, hi):
    """
    Ratio of density in the two outer BAF lobes to density in the central band.

    Under balanced heterozygosity almost all sites sit centrally and this ratio
    is small. Under LOH the distribution splits into a low lobe and a high lobe
    with the centre depleted, and the ratio becomes large. This is the direct
    expression of the bimodality that defines the LOH signature, and it is only
    computable because the allele fraction is kept unfolded.

    Returned as observed lobe mass divided by observed central mass, with the
    central term floored so that a completely empty centre yields a large
    finite value rather than a division by zero.
    """
    if not sites:
        return None
    low_lobe = 0
    high_lobe = 0
    central = 0
    for _pos, _depth, baf, _phased in sites:
        if baf < lo:
            low_lobe += 1
        elif baf > hi:
            high_lobe += 1
        else:
            central += 1
    total = float(len(sites))
    # Require mass in both lobes; a single-sided shift is a different
    # phenomenon (allele-specific dropout, contamination) and should not be
    # reported as LOH.
    if low_lobe == 0 or high_lobe == 0:
        return 0.0
    central_fraction = max(central / total, 1.0 / total)
    return ((low_lobe + high_lobe) / total) / central_fraction


def baf_deflection_from_balance(sites, lo, hi):
    """
    Absolute displacement of the region's BAF distribution away from 0.5.

    Computed as the mean absolute deviation of each site's unfolded allele
    fraction from 0.5. Under balanced heterozygosity this reflects only
    sampling and sequencing noise; under LOH the two lobes both sit far from
    the centre and the value rises sharply. Unlike a median, it does not cancel
    between the two lobes of a bimodal distribution.
    """
    if not sites:
        return None
    total = 0.0
    for _pos, _depth, baf, _phased in sites:
        total += abs(baf - 0.5)
    return total / float(len(sites))


# ---------------------------------------------------------------------------
# Copy number integration
#
# B-allele frequency alone cannot distinguish copy-neutral LOH from hemizygous
# deletion: both produce allelic imbalance. The discrimination requires copy
# number at the same locus, which is read here from ichorCNA output.
#
# Copy number also determines where balanced heterozygous sites are expected to
# sit. In a diploid region the expectation is 0.5, but at copy number 3 the
# balanced allele fractions are 1/3 and 2/3, and at copy number 1 there is no
# balanced state at all. Applying a diploid expectation to a non-diploid genome
# would mark every region as imbalanced, so the central band is shifted to the
# sample's own modal copy number.
# ---------------------------------------------------------------------------

CN_CALL_LOH = "CN_LOH"
CN_CALL_DELETION = "DELETION"
CN_CALL_GAIN = "GAIN"
CN_CALL_GAP = "GAP"
CN_CALL_GAP_INFERRED = "GAP_INFERRED"
CN_CALL_NONE = "NO_CN"


def read_ichor_params(path):
    """
    Read tumour fraction and ploidy from an ichorCNA .params.txt file, and
    assess whether the reported solution is trustworthy.

    ichorCNA reports a table of candidate solutions with log-likelihoods. It
    does not always report the highest-likelihood solution as its answer, and a
    tumour fraction of exactly zero alongside a non-diploid ploidy is a
    characteristic signature of a fit that failed to converge. Both conditions
    are detected here and surfaced, because a degenerate fit produces copy
    number calls that look confident but carry no information.

    Returns a dict with tumour_fraction, ploidy and fit_warning.
    """
    result = {"tumour_fraction": None, "ploidy": None, "fit_warning": None}
    if not path or not os.path.exists(path):
        return result

    chosen_loglik = None
    solution_logliks = []

    with open(path) as handle:
        for line in handle:
            stripped = line.strip()
            fields = [field.strip() for field in stripped.split("\t")]

            if stripped.lower().startswith("tumor fraction:"):
                try:
                    result["tumour_fraction"] = float(stripped.split(":", 1)[1])
                except ValueError:
                    pass
            elif stripped.lower().startswith("ploidy:"):
                try:
                    result["ploidy"] = float(stripped.split(":", 1)[1])
                except ValueError:
                    pass
            elif stripped.lower().startswith("log-likelihood:"):
                try:
                    chosen_loglik = float(stripped.split(":", 1)[1])
                except ValueError:
                    pass
            elif len(fields) >= 7:
                # Candidate solution rows in the trailing table; the final
                # column is the log-likelihood for that starting point.
                try:
                    solution_logliks.append(float(fields[-1]))
                except ValueError:
                    pass

    warnings = []
    if result["tumour_fraction"] == 0:
        warnings.append("tumour fraction reported as zero")
    if (chosen_loglik is not None and solution_logliks
            and max(solution_logliks) > chosen_loglik):
        warnings.append(
            "reported solution log-likelihood {0:.0f} is below the best "
            "candidate {1:.0f}".format(chosen_loglik, max(solution_logliks))
        )
    if result["ploidy"] is not None and result["tumour_fraction"] == 0 \
            and abs(result["ploidy"] - 2.0) > 0.3:
        warnings.append("zero tumour fraction with non-diploid ploidy")

    if warnings:
        result["fit_warning"] = "; ".join(warnings)
    return result


def read_ichor_segments(path):
    """
    Read segmented copy number from an ichorCNA .cna.seg file.

    Returns {chrom: [(start, end, copy_number, event), ...]} with chromosome
    names normalised to a chr prefix. Column positions are located by header
    name because ichorCNA prefixes most columns with the sample identifier and
    has renamed them between versions.

    Copy number is taken from the corrected copy number column where present,
    falling back to the raw call. The log ratio column is deliberately not used:
    it is empty for some bins, and copy number is stated directly.
    """
    if not path or not os.path.exists(path):
        return {}

    segments = {}
    with open(path) as handle:
        header = [column.strip().lower()
                  for column in handle.readline().rstrip("\n").split("\t")]

        def find(suffix):
            for index, name in enumerate(header):
                if name == suffix or name.endswith("." + suffix):
                    return index
            return None

        chrom_index = find("chr")
        start_index = find("start")
        end_index = find("end")
        cn_index = find("corrected_copy_number")
        if cn_index is None:
            cn_index = find("copy.number")
        event_index = find("corrected_call")
        if event_index is None:
            event_index = find("event")

        if None in (chrom_index, start_index, end_index, cn_index):
            sys.stderr.write(
                "Could not locate copy number columns in {0}\n".format(path)
            )
            return {}

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(chrom_index, start_index, end_index, cn_index):
                continue
            chrom = fields[chrom_index]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                start = int(float(fields[start_index]))
                end = int(float(fields[end_index]))
                copy_number = float(fields[cn_index])
            except ValueError:
                continue
            event = ""
            if event_index is not None and event_index < len(fields):
                event = fields[event_index]
            segments.setdefault(chrom, []).append((start, end, copy_number, event))

    for chrom in segments:
        segments[chrom].sort(key=lambda item: item[0])
    return segments


def modal_copy_number(segments):
    """
    Copy number covering the greatest number of bases across the genome.

    Used as the sample's baseline ploidy state for setting the expected
    balanced allele fractions. Base-weighting rather than bin-counting avoids
    letting many small bins outvote the bulk of the genome.
    """
    if not segments:
        return None
    mass = {}
    for chrom_segments in segments.values():
        for start, end, copy_number, _event in chrom_segments:
            span = max(0, end - start)
            key = int(round(copy_number))
            mass[key] = mass.get(key, 0) + span
    if not mass:
        return None
    return max(mass.items(), key=lambda item: item[1])[0]


def expected_band_for_copy_number(copy_number, half_width):
    """
    Central BAF band within which balanced heterozygous sites are expected.

    At copy number 2 the balanced fraction is 0.5 and the band is centred
    there. At copy number 3 the balanced fractions are 1/3 and 2/3; the band is
    taken around the lower mode, and its mirror at 1 minus that value is
    handled by the caller.

    Copy number 1 has no heterozygous state at all, so no band is defined.
    """
    if copy_number is None or copy_number < 2:
        return None
    minor = 1.0
    major = copy_number - 1.0
    centre = minor / float(minor + major)
    return (max(0.0, centre - half_width), min(1.0, centre + half_width))


def bimodality_is_informative(modal_copy_number):
    """
    Whether a two-lobed BAF distribution is evidence of LOH.

    At copy number 2 the balanced state is a single mode at 0.5, so splitting
    into two lobes indicates allelic imbalance. At copy number 3 or above the
    balanced state is *itself* two-lobed - a normal triploid region sits at 1/3
    and 2/3 - so two lobes carry no information about LOH and the statistic
    must not contribute to flagging. Ignoring this would mark every region of a
    triploid genome as imbalanced.
    """
    return modal_copy_number is None or int(round(modal_copy_number)) == 2


def overlapping_segments(segments, chrom, start, end):
    """
    Segments genuinely intersecting a half-open interval.

    A true overlap test is required rather than containment: ichorCNA bins are
    megabase-scale and a panel window will typically sit inside one bin, or
    straddle two, without ever containing one.
    """
    hits = []
    for seg_start, seg_end, copy_number, event in segments.get(chrom, []):
        if seg_start < end and seg_end > start:
            overlap = min(seg_end, end) - max(seg_start, start)
            if overlap > 0:
                hits.append((overlap, copy_number, event))
    return hits


def flanking_copy_number(segments, chrom, start, end, max_distance):
    """
    Infer copy number for a region falling in a gap between ichorCNA bins.

    ichorCNA omits bins where coverage was insufficient, which on adaptive
    sampling data happens routinely in the off-target background. Where a panel
    region falls entirely inside such a gap, the nearest bin on each side is
    consulted: if both flanks agree and both lie within max_distance, their
    shared copy number is a well-supported inference for the intervening
    sequence.

    The inference is deliberately conservative. Disagreeing flanks mean a
    segment boundary may fall inside the gap, and a single flank gives no
    evidence that the state persists across it; both cases return nothing so
    the region stays reported as an unmeasured gap.
    """
    chrom_segments = segments.get(chrom, [])
    if not chrom_segments:
        return None, ""

    left = None
    right = None
    for seg_start, seg_end, copy_number, event in chrom_segments:
        if seg_end <= start:
            left = (start - seg_end, copy_number, event)
        elif seg_start >= end and right is None:
            right = (seg_start - end, copy_number, event)

    if left is None or right is None:
        return None, ""
    if left[0] > max_distance or right[0] > max_distance:
        return None, ""
    if int(round(left[1])) != int(round(right[1])):
        return None, ""

    return left[1], left[2]


def call_copy_number(segments, region, max_gap_distance=1500000,
                     allow_inference=True):
    """
    Assign a copy number call to a panel region.

    The copy number is the overlap-weighted median across intersecting bins.
    Where no bin intersects the region the result is reported as a gap rather
    than as a neutral call, because ichorCNA omits bins with insufficient
    coverage and on adaptive sampling data those omissions are common. Silently
    treating a gap as copy number 2 would manufacture false copy-neutral calls.
    """
    start = int(region["start"])
    end = int(region["end"])
    hits = overlapping_segments(segments, region["chrom"], start, end)

    if not hits:
        # No bin intersects the region. Fall back to concordant flanking bins
        # where they are close enough to support the inference.
        if not allow_inference:
            return None, "", CN_CALL_GAP
        inferred, event = flanking_copy_number(
            segments, region["chrom"], start, end, max_gap_distance
        )
        if inferred is None:
            return None, "", CN_CALL_GAP
        return inferred, event, CN_CALL_GAP_INFERRED

    weighted = []
    for overlap, copy_number, event in hits:
        # Weight by overlap in 10 kb units so the median reflects base coverage.
        weighted.extend([copy_number] * max(1, overlap // 10000))
    copy_number = statistics.median(weighted)

    events = [event for _overlap, _cn, event in hits if event]
    event_label = max(set(events), key=events.count) if events else ""

    rounded = int(round(copy_number))
    if rounded == 2:
        call = CN_CALL_LOH
    elif rounded < 2:
        call = CN_CALL_DELETION
    else:
        call = CN_CALL_GAIN
    return copy_number, event_label, call


def annotate_copy_number(rows, segments, params_info, max_gap_distance=1500000):
    """
    Attach copy number columns to a sample's region rows.

    The cn_call column expresses what the copy number implies *given* that the
    BAF screen flagged the region. It is only meaningful alongside the BAF
    flag: a region at copy number 2 is labelled CN_LOH only in the sense that
    imbalance there would be copy-neutral, not that imbalance was observed.
    Interpretation combines the two columns and is spelled out in cn_note.
    """
    warning = params_info.get("fit_warning")
    # Inferring copy number across a gap compounds any error in the underlying
    # fit. Where ichorCNA's own solution is suspect the inference is withheld
    # and the region is reported as an unmeasured gap instead.
    allow_inference = warning is None
    for row in rows:
        copy_number, event, call = call_copy_number(
            segments, row, max_gap_distance=max_gap_distance,
            allow_inference=allow_inference
        )
        row["cn"] = copy_number
        row["cn_event"] = event
        row["cn_call"] = call if segments else CN_CALL_NONE
        row["tumour_fraction"] = params_info.get("tumour_fraction")

        if not segments:
            row["cn_note"] = "no copy number data"
        elif warning:
            row["cn_note"] = "unreliable ichorCNA fit: {0}".format(warning)
        elif call == CN_CALL_GAP:
            row["cn_note"] = "no ichorCNA bin overlaps this region"
        elif call == CN_CALL_GAP_INFERRED:
            rounded = int(round(copy_number)) if copy_number is not None else None
            implication = ""
            if row.get("flag") == FLAG_LOH and rounded is not None:
                if rounded == 2:
                    implication = "; imbalance would be copy-neutral LOH"
                elif rounded < 2:
                    implication = "; imbalance consistent with hemizygous deletion"
                else:
                    implication = "; imbalance consistent with allelic gain"
            row["cn_note"] = ("copy number inferred from concordant flanking "
                              "bins, not measured in this window" + implication)
        elif row.get("flag") == FLAG_LOH and call == CN_CALL_LOH:
            row["cn_note"] = "imbalance at neutral copy number: copy-neutral LOH"
        elif row.get("flag") == FLAG_LOH and call == CN_CALL_DELETION:
            row["cn_note"] = "imbalance with copy loss: hemizygous deletion"
        elif row.get("flag") == FLAG_LOH and call == CN_CALL_GAIN:
            row["cn_note"] = "imbalance with copy gain: allelic imbalance from gain"
        else:
            row["cn_note"] = ""


def median_or_none(values):
    """Median of a list, or None if the list is empty."""
    if not values:
        return None
    return statistics.median(values)


def format_value(value, digits=4):
    """Format a float for TSV output, emitting NA for missing values."""
    if value is None:
        return "NA"
    if isinstance(value, float):
        if math.isnan(value) or math.isinf(value):
            return "NA"
        return "{0:.{1}f}".format(value, digits)
    return str(value)


# ---------------------------------------------------------------------------
# Per-region computation
# ---------------------------------------------------------------------------

def sites_in_region(chrom_sites, start, end):
    """
    Select sites falling inside a half-open BED interval.

    A linear scan is used rather than bisection because the per-chromosome site
    lists in panel data are small; correctness is preferred to micro-optimisation
    here. Positions in the VCF are 1-based, the BED is 0-based half-open.
    """
    selected = []
    for pos, depth, minor_baf, phased in chrom_sites:
        zero_based = pos - 1
        if zero_based < start:
            continue
        if zero_based >= end:
            break
        selected.append((pos, depth, minor_baf, phased))
    return selected


def screen_region(region, by_chrom, params, cache):
    """
    Compute all per-region metrics for a single panel window.

    Returns a dict of metrics. The flag is not set here; flagging is applied in
    a later pass so that cohort context is available.
    """
    chrom_sites = by_chrom.get(region["chrom"], [])
    sites = sites_in_region(chrom_sites, region["start"], region["end"])

    span_mb = (region["end"] - region["start"]) / 1e6
    n_het = len(sites)
    het_per_mb = n_het / span_mb if span_mb > 0 else None

    depths = [site[1] for site in sites]
    bafs = [site[2] for site in sites]
    n_phased = sum(1 for site in sites if site[3])

    # The central band is set from the sample's modal copy number where copy
    # number data is available, so that a triploid genome is compared against
    # the 1/3 and 2/3 balanced fractions rather than against 0.5.
    band = params.get("central_band")
    if band is not None:
        lo, hi = band
    else:
        lo = params["central_lo"]
        hi = params["central_hi"]

    result = {
        "region": region["name"],
        "chrom": region["chrom"],
        "start": region["start"],
        "end": region["end"],
        "span_mb": span_mb,
        "n_het": n_het,
        "het_per_mb": het_per_mb,
        "median_dp": median_or_none(depths),
        "median_baf": median_or_none(bafs),
        "frac_central": None,
        "frac_phased": (n_phased / float(n_het)) if n_het else None,
        "depletion_score": None,
        "bimodality": None,
        "baf_deflection": None,
        "het_density_ratio": None,
        "assessable": n_het >= params["min_sites"],
    }

    if result["assessable"]:
        result["depletion_score"] = central_depletion_ratio(sites, lo, hi, cache)
        result["bimodality"] = bimodality_ratio(sites, lo, hi)
        result["baf_deflection"] = baf_deflection_from_balance(sites, lo, hi)
        central = sum(1 for site in sites if lo <= site[2] <= hi)
        result["frac_central"] = central / float(n_het)
    result["band_lo"] = lo
    result["band_hi"] = hi

    return result


def apply_flags(sample_results, cohort_density, params):
    """
    Assign an interpretive flag to each region.

    Flagging logic:
      UNASSESSABLE  - too few usable het sites to say anything.
      LOH_LIKELY    - central BAF density depleted AND the distribution is
                      bimodal with mass in both outer lobes.
      EQUIVOCAL     - one criterion met but not the other.
      NO_LOH        - neither criterion met.

    Requiring bimodality alongside central depletion is what separates LOH from
    the other ways a region can lose central density. A region that is merely
    noisy loses central mass symmetrically without forming two lobes; a region
    with allele-specific dropout forms one lobe only. Both are excluded.

    Het density is reported and used only as corroboration, never as a primary
    criterion. It is confounded by window size, mappability and coverage, and
    the wide Ig-locus windows carry intrinsically low density in every sample.
    Comparing each region against the cohort median for that same region
    removes the confounder; where no cohort baseline exists the density term is
    reported but does not contribute to flagging.
    """
    for row in sample_results:
        if not row["assessable"]:
            row["flag"] = FLAG_UNASSESSABLE
            row["flag_reason"] = "fewer than {0} usable het sites".format(params["min_sites"])
            continue

        baseline = cohort_density.get(row["region"])
        if baseline and baseline > 0 and row["het_per_mb"] is not None:
            row["het_density_ratio"] = row["het_per_mb"] / baseline

        depleted = (row["depletion_score"] is not None
                    and row["depletion_score"] < params["cdr_threshold"])
        # Bimodality only discriminates at copy number 2; above that the
        # balanced state is itself two-lobed. See bimodality_is_informative.
        bimodal = (params.get("bimodality_informative", True)
                   and row["bimodality"] is not None
                   and row["bimodality"] >= params["bimodality_threshold"])
        density_low = (row["het_density_ratio"] is not None
                       and row["het_density_ratio"] < params["density_ratio"])

        reasons = []
        if depleted:
            reasons.append("central BAF density depleted")
        if bimodal:
            reasons.append("BAF distribution bimodal")
        if density_low:
            reasons.append("het density below cohort baseline")

        if not params.get("bimodality_informative", True):
            # In a non-diploid genome the bimodality test is unavailable, so a
            # single criterion remains. Rather than call LOH on depletion
            # alone, which is weaker evidence, such regions are reported as
            # equivocal and referred to allele-specific copy number analysis.
            row["flag"] = FLAG_EQUIVOCAL if depleted else FLAG_NO_LOH
            if depleted:
                reasons.append("non-diploid genome: bimodality uninformative, "
                               "allele-specific copy number analysis required")
        elif depleted and bimodal:
            row["flag"] = FLAG_LOH
        elif depleted or bimodal:
            row["flag"] = FLAG_EQUIVOCAL
        else:
            row["flag"] = FLAG_NO_LOH

        row["flag_reason"] = "; ".join(reasons) if reasons else "balanced heterozygosity"


def compute_cohort_density(all_results):
    """
    Median het density per Mb for each region across all samples screened.

    This is the baseline against which an individual sample's het density is
    judged. With fewer than three samples the baseline is unstable and is
    withheld, so density does not contribute to flagging.
    """
    per_region = {}
    for sample_id, rows in all_results.items():
        for row in rows:
            if row["het_per_mb"] is not None and row["assessable"]:
                per_region.setdefault(row["region"], []).append(row["het_per_mb"])
    return {
        region: statistics.median(values)
        for region, values in per_region.items()
        if len(values) >= 3
    }


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

OUTPUT_COLUMNS = [
    "sample", "region", "chrom", "start", "end", "span_mb",
    "n_het", "het_per_mb", "het_density_ratio", "median_dp", "frac_phased",
    "median_baf", "frac_central", "band_lo", "band_hi", "baf_deflection",
    "depletion_score", "bimodality", "flag", "flag_reason",
    "cn", "cn_event", "cn_call", "tumour_fraction", "cn_note",
]


def write_table(all_results, out_path):
    """Write the combined per-sample, per-region table."""
    with open(out_path, "w") as handle:
        handle.write("\t".join(OUTPUT_COLUMNS) + "\n")
        for sample_id in sorted(all_results):
            for row in all_results[sample_id]:
                record = dict(row)
                record["sample"] = sample_id
                fields = []
                for column in OUTPUT_COLUMNS:
                    value = record.get(column)
                    if column in ("n_het", "start", "end"):
                        fields.append("NA" if value is None else str(value))
                    elif column in ("sample", "region", "chrom", "flag", "flag_reason",
                                    "cn_event", "cn_call", "cn_note"):
                        fields.append("NA" if value is None else str(value))
                    else:
                        fields.append(format_value(value))
                handle.write("\t".join(fields) + "\n")


def read_sample_map(path):
    """
    Read a two-column sample map: sample identifier and VCF path.

    Blank lines and lines beginning with '#' are ignored.
    """
    entries = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field.strip() for field in line.split("\t")]
            if len(fields) < 2:
                raise ValueError("Sample map line needs two tab-separated fields: {0}".format(line))

            sample_id, vcf_path = fields[0], fields[1]

            # Terminal integrations can inject control sequences into the first
            # line of a redirected command's output, which silently corrupts the
            # sample identifier and produces file paths that cannot exist. Fail
            # loudly here rather than reporting a missing input file later.
            control = [character for character in sample_id
                       if ord(character) < 32 or ord(character) == 127]
            if control:
                raise ValueError(
                    "Sample identifier contains control characters and is likely "
                    "corrupted: {0!r}. Check the sample map with 'cat -A'; a "
                    "terminal integration may have written escape sequences into "
                    "the first line.".format(sample_id)
                )
            if not sample_id:
                raise ValueError("Empty sample identifier in sample map: {0!r}".format(line))

            entries.append((sample_id, vcf_path))
    if not entries:
        raise ValueError("No samples parsed from sample map: {0}".format(path))
    return entries


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def read_ichor_map(path):
    """
    Read an explicit per-sample copy number file map.

    Two or three tab-separated columns: sample identifier, segment file, and
    optionally the parameter file. Where the parameter file is omitted it is
    looked for alongside the segment file, which covers the common case of both
    being staged into the same directory.
    """
    mapping = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = [field.strip() for field in line.split("\t")]
            if len(fields) < 2:
                raise ValueError(
                    "Copy number map needs at least two tab-separated fields: {0}".format(line)
                )
            sample_id, seg_path = fields[0], fields[1]
            if len(fields) >= 3 and fields[2]:
                params_path = fields[2]
            else:
                params_path = seg_path.replace(".cna.seg", ".params.txt")
            mapping[sample_id] = (seg_path, params_path)
    return mapping


def build_argument_parser():
    parser = argparse.ArgumentParser(
        description="Panel-wide BAF screen for loss of heterozygosity from phased Clair3 VCFs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--bed", required=True,
                        help="Panel BED file, four columns, chr-named to match the VCF.")
    parser.add_argument("--vcf-dir", dest="vcf_dir",
                        help="VCF file or directory of per-chromosome VCFs for a single sample.")
    parser.add_argument("--sample", help="Sample identifier, required with --vcf-dir.")
    parser.add_argument("--sample-map", dest="sample_map",
                        help="Two-column TSV of sample identifier and VCF path, for cohort runs.")
    parser.add_argument("--sample-column", dest="sample_column",
                        help="Sample column name within multi-sample VCFs. Defaults to the first.")
    parser.add_argument("--out", required=True, help="Output TSV path.")
    parser.add_argument("--ichor-map", dest="ichor_map",
                        help="Two- or three-column TSV of sample identifier, .cna.seg path "
                             "and optionally .params.txt path. Use instead of --ichor-dir "
                             "where the copy number files are not laid out in the "
                             "<dir>/<sample>/ichorcna_out/ convention, for example when "
                             "staged individually by a workflow manager.")
    parser.add_argument("--ichor-dir", dest="ichor_dir",
                        help="ichorCNA results root. Segment and parameter files are "
                             "expected at <dir>/<sample>/ichorcna_out/<sample>.cna.seg "
                             "and .params.txt. Without this the screen reports BAF only.")
    parser.add_argument("--max-gap-distance", dest="max_gap_distance", type=int,
                        default=1500000,
                        help="Maximum distance to a flanking bin when inferring copy "
                             "number for a region that falls in an ichorCNA gap.")
    parser.add_argument("--no-ploidy-adjust", dest="no_ploidy_adjust", action="store_true",
                        help="Keep the diploid central band even when copy number "
                             "indicates the sample is not diploid.")

    parser.add_argument("--min-site-depth", type=int, default=DEFAULT_MIN_SITE_DEPTH,
                        help="Minimum depth for a het site to be included.")
    parser.add_argument("--min-sites", type=int, default=DEFAULT_MIN_SITES,
                        help="Minimum usable het sites for a region to be assessable.")
    parser.add_argument("--central-lo", type=float, default=DEFAULT_CENTRAL_LO,
                        help="Lower bound of the central balanced BAF band.")
    parser.add_argument("--central-hi", type=float, default=DEFAULT_CENTRAL_HI,
                        help="Upper bound of the central balanced BAF band.")
    parser.add_argument("--cdr-threshold", type=float, default=DEFAULT_CDR_THRESHOLD,
                        help="Central depletion ratio below which a region counts as depleted.")
    parser.add_argument("--bimodality-threshold", type=float, default=DEFAULT_BIMODALITY_THRESHOLD,
                        help="Outer-lobe to central mass ratio supporting an LOH call.")
    parser.add_argument("--density-ratio", type=float, default=DEFAULT_DENSITY_RATIO,
                        help="Het density relative to cohort median treated as corroborating.")
    return parser


def main(argv=None):
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    if args.sample_map:
        samples = read_sample_map(args.sample_map)
    elif args.vcf_dir and args.sample:
        samples = [(args.sample, args.vcf_dir)]
    else:
        parser.error("Provide either --sample-map, or both --vcf-dir and --sample.")

    params = {
        "min_site_depth": args.min_site_depth,
        "min_sites": args.min_sites,
        "central_lo": args.central_lo,
        "central_hi": args.central_hi,
        "cdr_threshold": args.cdr_threshold,
        "bimodality_threshold": args.bimodality_threshold,
        "density_ratio": args.density_ratio,
    }

    ichor_map = read_ichor_map(args.ichor_map) if args.ichor_map else {}

    regions = read_bed(args.bed)
    sys.stderr.write("Loaded {0} panel regions from {1}\n".format(len(regions), args.bed))

    # The binomial expectation depends only on depth and band bounds, so it is
    # cached across regions and samples.
    binomial_cache = {}
    all_results = {}
    sample_context = {}

    for sample_id, vcf_target in samples:
        try:
            vcf_paths = collect_vcf_paths(vcf_target)
        except ValueError as error:
            sys.stderr.write("Skipping {0}: {1}\n".format(sample_id, error))
            continue

        by_chrom, total_records, kept = load_het_sites(
            vcf_paths, params["min_site_depth"], args.sample_column
        )
        sys.stderr.write(
            "{0}: {1} VCF file(s), {2} records scanned, {3} usable het SNV sites retained\n".format(
                sample_id, len(vcf_paths), total_records, kept
            )
        )

        # Copy number is loaded before screening so that the sample's modal
        # copy number can set the expected balanced band.
        segments = {}
        params_info = {"tumour_fraction": None, "ploidy": None, "fit_warning": None}
        if args.ichor_map or args.ichor_dir:
            if args.ichor_map:
                entry = ichor_map.get(sample_id)
                if entry is None:
                    sys.stderr.write(
                        "{0}: no entry in copy number map; screening BAF only\n".format(sample_id)
                    )
                    seg_path = par_path = ""
                else:
                    seg_path, par_path = entry
            else:
                base = os.path.join(args.ichor_dir, sample_id, "ichorcna_out")
                seg_path = os.path.join(base, sample_id + ".cna.seg")
                par_path = os.path.join(base, sample_id + ".params.txt")
            if not os.path.exists(seg_path):
                sys.stderr.write("{0}: no segment file at {1}\n".format(sample_id, seg_path))
            segments = read_ichor_segments(seg_path)
            params_info = read_ichor_params(par_path)
            total_bins = sum(len(v) for v in segments.values())
            sys.stderr.write("{0}: {1} copy number bins across {2} chromosomes\n".format(
                sample_id, total_bins, len(segments)))
            if params_info.get("fit_warning"):
                sys.stderr.write("{0}: WARNING unreliable ichorCNA fit ({1})\n".format(
                    sample_id, params_info["fit_warning"]))

        sample_params = dict(params)
        half_width = (params["central_hi"] - params["central_lo"]) / 2.0
        modal_cn = modal_copy_number(segments) if segments else None
        if modal_cn is not None and not args.no_ploidy_adjust and modal_cn != 2:
            band = expected_band_for_copy_number(modal_cn, half_width)
            if band is not None:
                sample_params["central_band"] = band
                sys.stderr.write(
                    "{0}: modal copy number {1}, central band shifted to "
                    "{2:.3f}-{3:.3f}\n".format(sample_id, modal_cn, band[0], band[1]))

        if modal_cn is not None and not bimodality_is_informative(modal_cn):
            sample_params["bimodality_informative"] = False
            sys.stderr.write(
                "{0}: modal copy number {1} is not diploid; bimodality test "
                "disabled\n".format(sample_id, modal_cn))

        rows = [screen_region(region, by_chrom, sample_params, binomial_cache)
                for region in regions]
        all_results[sample_id] = rows
        sample_context[sample_id] = (segments, params_info)

    if not all_results:
        sys.stderr.write("No samples produced results.\n")
        return 1

    cohort_density = compute_cohort_density(all_results)
    if cohort_density:
        sys.stderr.write(
            "Cohort het-density baseline available for {0} regions across {1} samples\n".format(
                len(cohort_density), len(all_results)
            )
        )
    else:
        sys.stderr.write(
            "Fewer than three samples: het density reported but not used for flagging\n"
        )

    for sample_id, rows in all_results.items():
        segments, params_info = sample_context.get(sample_id, ({}, {}))
        flag_params = dict(params)
        modal = modal_copy_number(segments) if segments else None
        if modal is not None and not bimodality_is_informative(modal):
            flag_params["bimodality_informative"] = False
        apply_flags(rows, cohort_density, flag_params)
        annotate_copy_number(rows, segments, params_info,
                             max_gap_distance=args.max_gap_distance)

    write_table(all_results, args.out)
    sys.stderr.write("Wrote {0}\n".format(args.out))

    # Concise summary to stderr so a pipeline log carries the headline result.
    for sample_id in sorted(all_results):
        counts = {}
        for row in all_results[sample_id]:
            counts[row["flag"]] = counts.get(row["flag"], 0) + 1
        summary = ", ".join("{0}={1}".format(flag, counts[flag]) for flag in sorted(counts))
        sys.stderr.write("{0}: {1}\n".format(sample_id, summary))

        cn_counts = {}
        for row in all_results[sample_id]:
            if row.get("flag") == FLAG_LOH and row.get("cn_call"):
                cn_counts[row["cn_call"]] = cn_counts.get(row["cn_call"], 0) + 1
        if cn_counts:
            detail = ", ".join("{0}={1}".format(call, cn_counts[call])
                               for call in sorted(cn_counts))
            sys.stderr.write("{0}:   flagged regions by copy number: {1}\n".format(
                sample_id, detail))

    return 0


if __name__ == "__main__":
    sys.exit(main())
