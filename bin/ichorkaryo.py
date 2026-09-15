#!/usr/bin/env python3
"""
ichorkaryo.py

Turn ichorCNA's segment output into the karyotype JSON the report already
consumes, so depth-based copy number comes from ichorCNA while the dashboard
contract stays unchanged.

Why this replaces mmkaryo's depth half
--------------------------------------
On a sample with FISH monosomy 13 at 95 percent and del(17p) at 95 percent,
mmkaryo returned a near-flat genome while ichorCNA called both:

    chr13  20-113 Mb   CN 1  logR -0.912  93 bins
    chr17   4-10  Mb   CN 1  logR -0.815   6 bins

The del(17p) is the instructive one. mmkaryo masks the adaptive sampling panel,
and the v7 panel window over TP53 is exactly the region the FISH probe targets,
so a focal 17p deletion is masked out of the data before segmentation. ichorCNA
bins the whole genome at 1 Mb without masking and sees it in six bins.

What is and is not carried over
-------------------------------
Available from ichorCNA: integer copy number, log2 ratio, segment bounds, the
HETD/NEUT/GAIN/AMP/HOMD call, and a subclonal flag — enough for arm calls, an
ISCN string, altered segments, and per-gene rows.

Not available: minor allele copy number and copy-neutral LOH. ichorCNA is
depth-only and cannot separate 2:0 from 1:1. That layer stays with mmbaf, whose
output is joined separately.

Purity
------
Taken from flow cytometry, not from ichorCNA. ichorCNA was built for cfDNA,
where tumour fraction means the proportion of circulating DNA that is tumour
and is genuinely small. These samples are sorted plasma cells at 92 to 99
percent, so there is no dilution to estimate: what ichorCNA reports is the
magnitude of copy-number deviation expressed as a fraction. On this cohort four
samples returned a tumour fraction near zero while carrying NRAS, ATM or MAX
variants at 35 to 69 percent allele fraction — copy-number-quiet genomes, not
tumour-poor ones. The ichorCNA value is carried through as cna_burden, named
for what it measures.

Arm calls
---------
ichorCNA segments do not respect arm boundaries — one sample has a single
segment spanning chr9:34-136 Mb across the centromere. Each arm therefore takes
the length-weighted majority call of the segments overlapping it, reported with
the fraction of the arm that call covers, so a partial event is visible as
partial rather than rounded into a whole-arm statement.

Cell fraction
-------------
Derived from the observed log2 against the expected log2 for a single-copy
change at the supplied purity. Assumes one clone; a segment ichorCNA flags
subclonal carries that flag alongside.

Dependencies: standard library only.
"""

import argparse
import csv
import json
import math
import os
import statistics
import sys
from collections import defaultdict

MAIN_CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

CALL_TO_EVENT = {
    "HOMD": "loss", "HETD": "loss",
    "NEUT": "none",
    "GAIN": "gain", "AMP": "gain", "HLAMP": "gain",
}


def read_cytobands(path):
    """Return {chrom: [(start, end, band, stain), ...]} for primary chromosomes."""
    bands = defaultdict(list)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom = parts[0]
            if chrom not in MAIN_CHROMS:
                continue
            bands[chrom].append((int(parts[1]), int(parts[2]), parts[3], parts[4]))
    for chrom in bands:
        bands[chrom].sort()
    return dict(bands)


def arm_bounds(bands):
    """Return {chrom: (p_start, p_end, q_start, q_end)} from the acen entries."""
    bounds = {}
    for chrom, entries in bands.items():
        acen = [(s, e) for s, e, _b, stain in entries if stain == "acen"]
        if not acen:
            continue
        cen_start = min(s for s, _e in acen)
        cen_end = max(e for _s, e in acen)
        chrom_start = min(s for s, _e, _b, _st in entries)
        chrom_end = max(e for _s, e, _b, _st in entries)
        bounds[chrom] = (chrom_start, cen_start, cen_end, chrom_end)
    return bounds


def band_for(bands, chrom, start, end):
    """Band label spanning a region, e.g. p13.2 or p13.2-p12."""
    if chrom not in bands:
        return ""
    hits = [b for s, e, b, _st in bands[chrom] if s < end and e > start]
    if not hits:
        return ""
    return hits[0] if len(hits) == 1 else f"{hits[0]}-{hits[-1]}"


def read_segments(path):
    """Read ichorCNA .seg.txt into a list of dicts."""
    segments = []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            chrom = row.get("chrom") or ""
            if chrom not in MAIN_CHROMS:
                continue
            try:
                start, end = int(row["start"]), int(row["end"])
                logr = float(row["seg.median.logR"])
            except (KeyError, ValueError):
                continue
            call = (row.get("Corrected_Call") or row.get("call") or "NEUT").strip()
            try:
                copies = int(float(row.get("Corrected_Copy_Number")
                                   or row.get("copy.number") or 2))
            except ValueError:
                copies = 2
            segments.append({
                "chrom": chrom, "start": start, "end": end,
                "logr": logr, "cn": copies, "call": call,
                "event": CALL_TO_EVENT.get(call, "none"),
                "n_bins": int(row.get("num.mark") or 0),
                "subclonal": str(row.get("subclone.status", "")).upper() == "TRUE",
            })
    return segments


def read_params(path):
    """Pull ploidy and ichorCNA's own fraction out of .params.txt."""
    out = {"ploidy": None, "cna_burden": None}
    if not path or not os.path.isfile(path):
        return out
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            parts = [p.strip() for p in line.rstrip("\n").split("\t")]
            if len(parts) < 2:
                continue
            key = parts[0].lower().rstrip(":")
            try:
                value = float(parts[1])
            except ValueError:
                continue
            if key.startswith("ploidy"):
                out["ploidy"] = value
            elif key.startswith("tumor fraction") or key.startswith("tumour fraction"):
                out["cna_burden"] = value
    return out


def expected_log2(copies, purity, expected_copies=2):
    """log2 ratio expected for an integer copy number at a given purity."""
    numerator = purity * copies + (1.0 - purity) * expected_copies
    denominator = 2.0
    if numerator <= 0:
        return float("-inf")
    return math.log2(numerator / denominator)


def cell_fraction(logr, copies, purity, expected_copies=2):
    """
    Fraction of cells carrying a single-copy change, from the observed log2.

    Observed copies are 2 * 2**logr; the difference from the expected copy
    number is the fraction of cells carrying the change, assuming one clone and
    a change of one copy. NaN outside (0, 1], which is where the change is more
    than one copy or the segment is not an event.
    """
    if copies == expected_copies:
        return None
    observed = 2.0 * (2.0 ** logr)
    value = abs(observed - expected_copies)
    if not math.isfinite(value) or value <= 0.0 or value > 1.0:
        return None
    return round(value, 3)


def expected_copies_for(chrom, sex, baseline=2):
    """Expected copy number for event calling, or None where it does not apply.

    Autosomes are judged against the genome's own baseline, not against 2. A
    tumour fitted at ploidy 3 whose arms sit uniformly at copy number 3 has a
    triploid baseline and no gains; comparing it to 2 reports every arm as a
    gain.

    The sex chromosomes keep a constitutional expectation rather than the
    tumour baseline, since their normal copy number is set by sex.
    """
    if chrom == "chrX":
        return 1 if sex == "XY" else 2
    if chrom == "chrY":
        return 1 if sex == "XY" else None
    return baseline


def noise_threshold(segments, baseline, floor=0.10, n_sigma=3.0):
    """Smallest |log2| worth calling: the larger of a floor and the noise.

    The noise term is n_sigma times the standard deviation of the log2 ratios
    of segments sitting at the baseline copy number, weighted by the number of
    bins each carries. Those segments should be flat, so their scatter is the
    measurement error for this sample.

    Falls back to the floor when there are too few baseline segments to
    estimate a spread, which is the case on a heavily rearranged genome.
    """
    values, weights = [], []
    for segment in segments:
        if segment["chrom"] in ("chrX", "chrY"):
            continue
        if segment["cn"] != baseline:
            continue
        bins = max(1, segment["n_bins"])
        values.append(segment["logr"])
        weights.append(bins)

    if len(values) < 3:
        return floor, None

    total = float(sum(weights))
    mean = sum(v * w for v, w in zip(values, weights)) / total
    variance = sum(w * (v - mean) ** 2 for v, w in zip(values, weights)) / total
    sigma = math.sqrt(variance) if variance > 0 else 0.0
    return max(floor, n_sigma * sigma), sigma


def autosomal_baseline(segments, bounds):
    """Modal autosomal copy number, weighted by the span each call covers.

    Length weighting matters: a genome can carry more distinct segments at one
    copy number while most of its sequence sits at another, and it is the
    sequence that defines the baseline.
    """
    covered = {}
    for segment in segments:
        chrom = segment["chrom"]
        if chrom in ("chrX", "chrY") or chrom not in bounds:
            continue
        span = max(0, segment["end"] - segment["start"])
        covered[segment["cn"]] = covered.get(segment["cn"], 0) + span
    if not covered:
        return 2, 0.0
    baseline = max(covered, key=lambda k: covered[k])
    total = sum(covered.values())
    return baseline, (covered[baseline] / total if total else 0.0)


def summarise_arms(segments, bounds, sex, baseline=2, threshold=0.0):
    """Length-weighted majority call per arm, with the fraction it covers."""
    arms = {}
    for chrom, (chrom_start, cen_start, cen_end, chrom_end) in sorted(bounds.items()):
        for label, arm_start, arm_end in (
            (chrom[3:] + "p", chrom_start, cen_start),
            (chrom[3:] + "q", cen_end, chrom_end),
        ):
            arm_len = arm_end - arm_start
            if arm_len <= 0:
                continue

            covered = defaultdict(int)
            weighted_logr = defaultdict(float)
            for segment in segments:
                if segment["chrom"] != chrom:
                    continue
                overlap = min(segment["end"], arm_end) - max(segment["start"], arm_start)
                if overlap <= 0:
                    continue
                covered[segment["cn"]] += overlap
                weighted_logr[segment["cn"]] += segment["logr"] * overlap

            if not covered:
                continue

            cn = max(covered, key=lambda k: covered[k])
            fraction = covered[cn] / arm_len
            logr = weighted_logr[cn] / covered[cn]
            expected = expected_copies_for(chrom, sex, baseline)

            if expected is None:
                event, fraction_cells = "none", None
            elif cn == expected:
                event, fraction_cells = "none", None
            elif abs(logr) < threshold:
                # The HMM assigned a copy state, but the depth does not differ
                # from baseline by enough to support it.
                event, fraction_cells = "none", None
            else:
                event = "gain" if cn > expected else "loss"
                fraction_cells = None  # filled below when purity is known

            arms[label] = {
                "chrom": chrom, "cn": cn, "log2": round(logr, 4),
                "event": event, "arm_fraction": round(fraction, 3),
                "cell_fraction": fraction_cells,
            }
    return arms


def iscn_string(arms, sex):
    """Group arms by copy number into an ISCN-style seq line."""
    by_cn = defaultdict(list)
    for label, value in arms.items():
        by_cn[value["cn"]].append(label)

    def sort_key(label):
        body = label.rstrip("pq")
        suffix = label[len(body):]
        index = {"X": 23, "Y": 24}.get(body)
        if index is None:
            try:
                index = int(body)
            except ValueError:
                index = 99
        return (index, {"p": 0, "q": 1, "": 2}.get(suffix, 3))

    parts = []
    for cn in sorted(by_cn):
        labels = sorted(by_cn[cn], key=sort_key)
        parts.append("(" + ", ".join(labels) + ")x%d" % cn)
    return "seq" + " ".join(parts)


def load_bed(path):
    """Read a panel BED into [(chrom, start, end, name), ...]."""
    regions = []
    if not path or not os.path.isfile(path):
        return regions
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            if len(parts) < 3 or parts[0] not in MAIN_CHROMS:
                continue
            name = parts[3] if len(parts) > 3 else ""
            regions.append((parts[0], int(parts[1]), int(parts[2]), name))
    return regions


CANONICAL_TRISOMY_CHROMS = ("chr3", "chr5", "chr7", "chr9",
                            "chr11", "chr15", "chr19", "chr21")


def read_cna_bins(path):
    """Per-bin copy number and log2 from ichorCNA's .cna.seg.

    Column names carry the sample id as a prefix, for example
    11F20265231.logR, so they are matched by suffix rather than exact name.
    Bins whose logR is NA come back with logr None: about 3 percent on this
    cohort, being the masked panel and centromere regions the run excludes by
    design.
    """
    bins = []
    if not path or not os.path.isfile(path):
        return bins
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []

        def column(suffix):
            for name in fields:
                if name == suffix or name.endswith("." + suffix):
                    return name
            return None

        logr_col = column("logR")
        cn_col = column("Corrected_Copy_Number") or column("copy.number")
        for row in reader:
            chrom = row.get("chr") or row.get("chrom") or ""
            if chrom not in MAIN_CHROMS:
                continue
            try:
                start, end = int(row["start"]), int(row["end"])
            except (KeyError, TypeError, ValueError):
                continue
            logr = None
            if logr_col:
                try:
                    logr = float(row[logr_col])
                except (TypeError, ValueError):
                    logr = None
                if logr is not None and logr != logr:
                    logr = None
            cn = None
            if cn_col:
                try:
                    cn = int(float(row[cn_col]))
                except (TypeError, ValueError):
                    cn = None
            bins.append({"chrom": chrom, "start": start, "end": end,
                         "logr": logr, "cn": cn})
    return bins


def bin_noise(bins, baseline, min_bins=30):
    """Per-bin log2 scatter at the baseline copy number, and the bin size.

    Deliberately not the segment-level sigma already computed for the event
    threshold. That one is the scatter of segment medians, smaller by roughly
    the square root of the segment length in bins, and it is the wrong
    quantity for asking whether a window of a given size could have been
    detected at all.

    Returns (sigma, bin_size, n_used). sigma is None when too few baseline
    bins carry a usable log2.
    """
    values, spans = [], []
    for entry in bins:
        # .cna.seg coordinates are 1-based inclusive, so a 1 Mb bin runs
        # 1000001 to 2000000 and its span is the difference plus one.
        span = entry["end"] - entry["start"] + 1
        if span > 1:
            spans.append(span)
        if entry["chrom"] in ("chrX", "chrY"):
            continue
        if entry["cn"] != baseline or entry["logr"] is None:
            continue
        values.append(entry["logr"])

    bin_size = int(statistics.median(spans)) if spans else None
    if len(values) < min_bins:
        return None, bin_size, len(values)

    # Median absolute deviation rather than the standard deviation: a handful
    # of bins at extreme corrected values would otherwise set the noise for
    # the whole genome.
    centre = statistics.median(values)
    mad = statistics.median([abs(v - centre) for v in values])
    sigma = 1.4826 * mad if mad > 0 else statistics.pstdev(values)
    return (sigma or None), bin_size, len(values)


def detection_limit(sigma, bin_size, span_bp, baseline, n_sigma=3.0):
    """Smallest single-copy change detectable over span_bp, as a cell fraction.

    A one-copy change present in a fraction f of cells shifts the mean log2
    over the window to log2((B - f) / B) on a baseline of B copies. It is
    detectable once that shift exceeds n_sigma standard errors of the window
    mean, which is sigma over the square root of the number of bins in the
    window.

    The baseline enters directly: on a triploid genome a single-copy loss is
    three to two, a smaller log2 step than two to one, and therefore harder to
    see. Returns None where the result is not a fraction of cells, which is
    the case for a window too small to resolve one copy at any fraction.
    """
    if not sigma or not bin_size or span_bp <= 0 or not baseline or baseline <= 0:
        return None
    n_bins = max(1.0, span_bp / float(bin_size))
    shift = n_sigma * sigma / math.sqrt(n_bins)
    fraction = baseline * (1.0 - 2.0 ** (-shift))
    if not math.isfinite(fraction) or fraction <= 0 or fraction > 1:
        return None
    return round(fraction, 4)


def read_wig(path):
    """Counts from a fixedStep WIG, as written by readCounter.

    Returns ({(chrom, start): count}, step), keyed by coordinate rather than
    pooled into a flat list, so the counts can be restricted to bins of a
    known copy number later. Coordinates follow the fixedStep header, which
    readCounter writes as start=1 with a fixed step, matching the start values
    in the .cna.seg.
    """
    counts, step, chrom, position = {}, None, None, None
    if not path or not os.path.isfile(path):
        return counts, step
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(("fixedStep", "variableStep", "track")):
                fields = dict(f.split("=", 1) for f in line.split()[1:]
                              if "=" in f)
                chrom = fields.get("chrom")
                try:
                    step = int(fields.get("step", step or 0)) or step
                    position = int(fields.get("start", 1))
                except ValueError:
                    position = 1
                continue
            if chrom is None or position is None or not step:
                continue
            try:
                counts[(chrom, position)] = float(line.split()[-1])
            except (IndexError, ValueError):
                pass
            position += step
    return counts, step


def wig_stats(counts_by_bin, cna_bins, baseline, min_bins=30):
    """Median reads per bin, and a robust overdispersion against Poisson.

    Overdispersion is measured only on bins the segmentation puts at the
    baseline copy number. Pooling the whole genome does not measure counting
    noise: a chromosome at three copies carries half again the reads of one at
    two, so the spread across a genome with any aberration is dominated by the
    aberration itself. Measured that way a synthetic genome with two trisomies
    and one monosomy returned 667, where the quantity is meant to sit near 1.

    Under counting noise alone the variance equals the mean and the ratio is
    1; the excess is everything else, which on adaptive sampling data is
    mostly mappability and GC. The variance comes from the median absolute
    deviation rather than the sample variance, since a few extreme bins would
    otherwise dominate.

    Returns (median_reads, overdispersion, n_baseline_bins). Overdispersion is
    None without a .cna.seg to say which bins are at baseline -- a wrong
    number here is worse than a missing one, since it is read as a data
    quality verdict.
    """
    positive = [c for c in counts_by_bin.values() if c > 0]
    if len(positive) < min_bins:
        return None, None, 0
    median_reads = int(round(statistics.median(positive)))

    at_baseline = []
    for entry in cna_bins:
        if entry["chrom"] in ("chrX", "chrY"):
            continue
        if entry["cn"] != baseline or entry["logr"] is None:
            continue
        count = counts_by_bin.get((entry["chrom"], entry["start"]))
        if count is not None and count > 0:
            at_baseline.append(count)

    if len(at_baseline) < min_bins:
        return median_reads, None, len(at_baseline)

    centre = statistics.median(at_baseline)
    if centre <= 0:
        return median_reads, None, len(at_baseline)
    mad = statistics.median([abs(v - centre) for v in at_baseline])
    sigma = 1.4826 * mad
    return median_reads, round((sigma * sigma) / centre, 2), len(at_baseline)


def chromosome_copy_numbers(segments, bounds):
    """Length-weighted majority copy number per whole chromosome."""
    out = {}
    for chrom, (chrom_start, _cen_start, _cen_end, chrom_end) in bounds.items():
        length = chrom_end - chrom_start
        if length <= 0:
            continue
        covered = defaultdict(int)
        for segment in segments:
            if segment["chrom"] != chrom:
                continue
            overlap = min(segment["end"], chrom_end) - max(segment["start"], chrom_start)
            if overlap > 0:
                covered[segment["cn"]] += overlap
        if not covered:
            continue
        cn = max(covered, key=lambda k: covered[k])
        out[chrom] = (cn, round(covered[cn] / length, 3))
    return out


def hyperdiploidy(chrom_cn, baseline=2, baseline_fraction=None,
                  min_fraction=0.8):
    """Trisomies and a ploidy class, against a constitutional two copies.

    The trisomy list is not referenced to the genome's modal baseline.
    Hyperdiploidy is defined against a normal diploid karyotype, and a sample
    whose modal copy number is three would otherwise report no trisomies at
    all. Most of a chromosome must sit at the trisomic level, so a large
    segmental gain is not read as a whole-chromosome trisomy.

    The hyperdiploid call does depend on the baseline, and that is the point.
    Two or more trisomies among 3, 5, 7, 9, 11, 15, 19 and 21 is the myeloma
    convention, but it only means hyperdiploidy against a diploid background.
    A uniformly near-triploid genome satisfies it trivially while being
    something else entirely: 11F20265231 returned fourteen autosomes at three
    copies with 83 percent of its span there, which is not a selective gain of
    odd chromosomes.

    near_triploid is kept separate from non_hyperdiploid rather than merged
    into it, since near-triploid myeloma is frequently a duplicated
    hyperdiploid clone and calling it non-hyperdiploid would be wrong in the
    other direction. The trisomy list is reported in every class.
    """
    def order(label):
        try:
            return int(label)
        except ValueError:
            return {"X": 23, "Y": 24}.get(label, 99)

    trisomies, autosome_copies, sex_copies = [], 0, 0
    for chrom, (cn, fraction) in chrom_cn.items():
        if chrom in ("chrX", "chrY"):
            sex_copies += cn
            continue
        autosome_copies += cn
        if cn == 3 and fraction >= min_fraction:
            trisomies.append(chrom[3:])

    canonical = [c for c in trisomies if "chr" + c in CANONICAL_TRISOMY_CHROMS]

    # An unreliable baseline makes every class below unreliable with it, so it
    # is named rather than guessed at.
    if (baseline_fraction is not None
            and baseline_fraction == baseline_fraction
            and baseline_fraction < 0.5):
        ploidy_class = "indeterminate"
    elif baseline <= 1:
        ploidy_class = "hypodiploid"
    elif baseline == 2:
        ploidy_class = ("hyperdiploid" if len(canonical) >= 2
                        else "non_hyperdiploid")
    elif baseline == 3:
        ploidy_class = "near_triploid"
    else:
        ploidy_class = "near_tetraploid"

    return {
        "ploidy_class": ploidy_class,
        "hyperdiploid": ploidy_class == "hyperdiploid",
        "trisomies": sorted(trisomies, key=order),
        "canonical_trisomies": sorted(canonical, key=order),
        "n_canonical_trisomies": len(canonical),
        "estimated_chromosome_count": autosome_copies + sex_copies,
    }


def corrected_overdispersion(sigma_log2, median_reads):
    """Per-bin log2 scatter over the scatter Poisson counting alone would give.

    The raw count overdispersion is taken before ichorCNA's GC and mappability
    correction and is dominated by both. This uses the corrected per-bin log2
    already computed for the detection limit, and compares it against the
    spread expected from counting N reads per bin, which for Poisson counts is
    1 / (ln2 * sqrt(N)) in log2 units.

    Near 1 means the measurement is counting-limited and more sequencing would
    help. Well above 1 means something other than read depth sets the noise
    floor, and more sequencing will not move it.
    """
    if not sigma_log2 or not median_reads or median_reads <= 0:
        return None
    expected = 1.0 / (math.log(2.0) * math.sqrt(median_reads))
    if expected <= 0:
        return None
    return round((sigma_log2 / expected) ** 2, 2)


def cytoband_table(bands, segments, purity, sex, baseline, threshold,
                   sigma, bin_size, min_fraction=0.5):
    """One row per altered cytoband, in the vocabulary a cytogenetics report uses.

    Bands are reported only where a single copy state covers at least
    min_fraction of the band, so a band straddling a breakpoint is left out
    rather than assigned to whichever side happens to be larger.
    """
    rows = []
    for chrom in sorted(bands, key=lambda c: MAIN_CHROMS.index(c)):
        expected = expected_copies_for(chrom, sex, baseline)
        if expected is None:
            continue
        for start, end, band, stain in bands[chrom]:
            span = end - start
            if span <= 0 or stain in ("acen", "gvar", "stalk"):
                continue
            covered = defaultdict(int)
            weighted = defaultdict(float)
            for segment in segments:
                if segment["chrom"] != chrom:
                    continue
                overlap = min(segment["end"], end) - max(segment["start"], start)
                if overlap <= 0:
                    continue
                covered[segment["cn"]] += overlap
                weighted[segment["cn"]] += segment["logr"] * overlap
            if not covered:
                continue
            cn = max(covered, key=lambda k: covered[k])
            fraction = covered[cn] / span
            if fraction < min_fraction or cn == expected:
                continue
            logr = weighted[cn] / covered[cn]
            if abs(logr) < threshold:
                continue
            rows.append({
                "cytoband": chrom[3:] + band,
                "chromosome": chrom,
                "start": start,
                "end": end,
                "copy_number": cn,
                "implied_copies": round(2.0 * (2.0 ** logr), 2),
                "log2": round(logr, 4),
                "event": "gain" if cn > expected else "loss",
                "band_fraction": round(fraction, 3),
                "cell_fraction": (cell_fraction(logr, cn, purity, expected)
                                  if purity else None),
                "detection_limit": detection_limit(sigma, bin_size, span, baseline),
            })
    return rows


def genes_from_panel(regions, segments, bands, purity, sex, baseline=2,
                     threshold=0.0):
    """One row per panel region, from the segment covering its midpoint."""
    rows = []
    for chrom, start, end, name in regions:
        if not name:
            continue
        midpoint = (start + end) // 2
        covering = [s for s in segments
                    if s["chrom"] == chrom and s["start"] <= midpoint < s["end"]]
        spanning = [s for s in segments
                    if s["chrom"] == chrom and s["start"] < end and s["end"] > start]

        if not covering:
            rows.append({"gene": name, "chrom": chrom, "start": start, "end": end,
                         "band": band_for(bands, chrom, start, end),
                         "cn": None, "log2": None, "event": "no_data",
                         "cell_fraction": None, "subclonal": False,
                         "split_across_segments": len(spanning) > 1})
            continue

        segment = covering[0]
        expected = expected_copies_for(chrom, sex, baseline)
        if expected is None or segment["cn"] == expected:
            event, fraction = "none", None
        elif abs(segment["logr"]) < threshold:
            event, fraction = "none", None
        else:
            event = "gain" if segment["cn"] > expected else "loss"
            fraction = cell_fraction(segment["logr"], segment["cn"], purity, expected) \
                if purity else None

        rows.append({
            "gene": name, "chrom": chrom, "start": start, "end": end,
            "band": band_for(bands, chrom, start, end),
            "cn": segment["cn"], "log2": round(segment["logr"], 4),
            "event": event, "cell_fraction": fraction,
            "subclonal": segment["subclonal"],
            "split_across_segments": len(spanning) > 1,
        })
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--seg", required=True, help="ichorCNA <sample>.seg.txt")
    parser.add_argument("--params", help="ichorCNA <sample>.params.txt")
    parser.add_argument("--cytobands", required=True, help="hg38 cytoBand file")
    parser.add_argument("--panel-bed", help="hg38 panel BED for this sample")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--purity", type=float,
                        help="flow cytometry purity as a fraction; without it no "
                             "cell fractions are reported")
    parser.add_argument("--sex", choices=["XX", "XY"], default=None)
    parser.add_argument("--min-log2", type=float, default=0.10,
                        help="Floor on the absolute log2 ratio required to call "
                             "an event. The effective threshold is the larger of "
                             "this and the sample's own noise. [default: 0.10]")
    parser.add_argument("--n-sigma", type=float, default=3.0,
                        help="Multiples of the baseline-segment standard "
                             "deviation used as the noise threshold. [default: 3]")
    parser.add_argument("--baseline", type=int, default=None,
                        help="Override the autosomal baseline copy number. By "
                             "default it is the modal arm copy number weighted "
                             "by span, which is right for a uniformly triploid "
                             "genome and wrong where no single level dominates.")
    parser.add_argument("--cna", help="ichorCNA <sample>.cna.seg. Supplies the "
                                      "per-bin log2 scatter the detection "
                                      "limits are derived from; without it "
                                      "they are null.")
    parser.add_argument("--wig", help="readCounter <sample>.wig. Supplies "
                                      "median reads per bin and "
                                      "overdispersion; without it both are "
                                      "null, since ichorCNA's own outputs "
                                      "carry log ratios and not counts.")
    parser.add_argument("--out", required=True, help="output karyotype JSON")
    args = parser.parse_args()

    if args.purity is not None and not 0.0 < args.purity <= 1.0:
        sys.exit(f"ERROR: purity {args.purity} is not a fraction in (0, 1]")

    bands = read_cytobands(args.cytobands)
    bounds = arm_bounds(bands)
    segments = read_segments(args.seg)
    params = read_params(args.params)

    if not segments:
        sys.exit(f"ERROR: no usable segments in {args.seg}")

    sex = args.sex or "XX"

    baseline, baseline_fraction = autosomal_baseline(segments, bounds)
    if args.baseline is not None:
        baseline = args.baseline
        baseline_fraction = float("nan")

    threshold, noise_sigma = noise_threshold(
        segments, baseline, floor=args.min_log2, n_sigma=args.n_sigma)

    arms = summarise_arms(segments, bounds, sex, baseline, threshold)

    if args.purity:
        for value in arms.values():
            if value["event"] != "none":
                expected = expected_copies_for(value["chrom"], sex, baseline)
                value["cell_fraction"] = cell_fraction(
                    value["log2"], value["cn"], args.purity, expected or baseline)

    altered = []
    n_below = 0
    for segment in segments:
        expected = expected_copies_for(segment["chrom"], sex, baseline)
        if expected is None or segment["cn"] == expected:
            continue
        if abs(segment["logr"]) < threshold:
            n_below += 1
            continue
        event = "gain" if segment["cn"] > expected else "loss"
        altered.append({
            "region": "%s:%d-%d" % (segment["chrom"], segment["start"], segment["end"]),
            "band": band_for(bands, segment["chrom"], segment["start"], segment["end"]),
            "size_mb": round((segment["end"] - segment["start"]) / 1e6, 2),
            "log2": round(segment["logr"], 4),
            "cn": segment["cn"],
            "call": segment["call"],
            "event": event,
            "n_bins": segment["n_bins"],
            "subclonal": segment["subclonal"],
            "cell_fraction": cell_fraction(segment["logr"], segment["cn"],
                                           args.purity, expected) if args.purity else None,
        })

    genes = genes_from_panel(load_bed(args.panel_bed), segments, bands,
                             args.purity, sex, baseline,
                             threshold) if args.panel_bed else []

    cna_bins = read_cna_bins(args.cna)
    sigma_bin, bin_size, n_baseline_bins = bin_noise(cna_bins, baseline)
    wig_counts, _wig_step = read_wig(args.wig)
    median_reads, overdispersion, n_wig_baseline = wig_stats(
        wig_counts, cna_bins, baseline)

    chrom_cn = chromosome_copy_numbers(segments, bounds)
    ploidy_calls = hyperdiploidy(chrom_cn, baseline, baseline_fraction)
    overdispersion_corrected = corrected_overdispersion(sigma_bin, median_reads)

    cytobands = cytoband_table(bands, segments, args.purity, sex, baseline,
                               threshold, sigma_bin, bin_size)

    # One list the report can render directly, rather than a single field the
    # parser has to know to look for.
    warnings = []
    if baseline_fraction == baseline_fraction and baseline_fraction < 0.5:
        warnings.append(
            "Less than half the autosomal span sits at the baseline copy "
            "number, so the modal value may not be the baseline. These calls "
            "need checking against the genome-wide plot.")
    if baseline != 2:
        warnings.append(
            "The modal copy number is %d, not 2. Gains and losses are called "
            "against %d, while trisomies and the hyperdiploid call are "
            "referenced to a constitutional 2 -- the two are not on the same "
            "scale and should not be read together." % (baseline, baseline))
    if ploidy_calls["ploidy_class"] == "near_triploid":
        warnings.append(
            "The genome is near-triploid: %d autosomes sit at three copies "
            "against a modal copy number of 3, so these are the baseline "
            "rather than gains. This is not hyperdiploidy in the myeloma "
            "sense, though a near-triploid clone is often a duplicated "
            "hyperdiploid one." % len(ploidy_calls["trisomies"]))
    elif ploidy_calls["ploidy_class"] == "near_tetraploid":
        warnings.append(
            "The modal copy number is %d. Gains, losses and the trisomy list "
            "should all be read against that rather than against a diploid "
            "genome." % baseline)
    if sigma_bin is None:
        warnings.append(
            "No per-bin log2 scatter available (%d usable baseline bins), so "
            "detection limits are not reported and a negative cannot be "
            "interpreted." % n_baseline_bins)
    if median_reads is None:
        warnings.append(
            "No read-count WIG supplied, so reads per bin and overdispersion "
            "are not reported.")
    elif overdispersion is None:
        warnings.append(
            "Reads per bin are reported but overdispersion is not: only %d "
            "bins could be matched to a baseline copy number, and measuring "
            "it across copy states would report the aberrations rather than "
            "the counting noise." % n_wig_baseline)

    summary = {
        "sample": args.sample,
        "source": "ichorCNA",
        "reference": "hg38",
        "fit": {
            "purity": args.purity,
            "purity_source": "flow" if args.purity else "not supplied",
            "ploidy": params["ploidy"],
            # Named for what it measures. ichorCNA's "tumour fraction" is the
            # magnitude of copy-number deviation; on sorted plasma cells it is
            # not the proportion of tumour present.
            "cna_burden": params["cna_burden"],
            "baseline_copy_number": baseline,
            "baseline_span_fraction": (None if baseline_fraction != baseline_fraction
                                       else round(baseline_fraction, 3)),
            "baseline_source": "supplied" if args.baseline is not None else "modal",
            "min_log2_threshold": round(threshold, 4),
            "noise_sigma_log2": (None if noise_sigma is None else round(noise_sigma, 4)),
            "segments_below_threshold": n_below,
            "baseline_warning": (
                "less than half the autosomal span sits at the baseline copy "
                "number; the modal value may not be the baseline and these "
                "calls need checking against the genome-wide plot"
                if baseline_fraction == baseline_fraction and baseline_fraction < 0.5
                else None),
        },
        "sex": sex,
        "n_segments": len(segments),
        "arms": arms,
        "ISCN_karyotype": iscn_string(arms, sex),
        "ploidy_class": ploidy_calls["ploidy_class"],
        "hyperdiploid": ploidy_calls["hyperdiploid"],
        "trisomies": ploidy_calls["trisomies"],
        "canonical_trisomies": ploidy_calls["canonical_trisomies"],
        "n_canonical_trisomies": ploidy_calls["n_canonical_trisomies"],
        "estimated_chromosome_count": ploidy_calls["estimated_chromosome_count"],
        "chromosome_copy_number": {c: {"cn": v[0], "fraction": v[1]}
                                   for c, v in sorted(chrom_cn.items())},
        "qc": {
            # Bin size here is readCounter's 1 Mb, not mmkaryo's 100 kb, so
            # reads per bin is not continuous with the previous numbers.
            "bin_size": bin_size,
            "n_bins": len(cna_bins),
            "n_bins_usable": sum(1 for b in cna_bins if b["logr"] is not None),
            "n_baseline_bins": n_baseline_bins,
            "sigma_log2_per_bin": (None if sigma_bin is None
                                   else round(sigma_bin, 4)),
            "detection_limit_5mb": detection_limit(sigma_bin, bin_size,
                                                   5_000_000, baseline),
            "median_reads_per_bin": median_reads,
            # Two different quantities, both on baseline bins only. The
            # corrected one is the one near 1 when counting-limited; the raw
            # one is before ichorCNA's GC and mappability correction and is
            # not comparable to mmkaryo's figure.
            "overdispersion": overdispersion_corrected,
            "overdispersion_corrected_log2": overdispersion_corrected,
            "overdispersion_raw_counts": overdispersion,
            "n_wig_baseline_bins": n_wig_baseline,
        },
        "cytobands": cytobands,
        "warnings": warnings,
        "altered_segments": altered,
        "genes": genes,
    }

    with open(args.out, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    n_events = sum(1 for v in arms.values() if v["event"] != "none")
    print(f"sample            : {args.sample}")
    print(f"segments          : {len(segments)}")
    print(f"altered segments  : {len(altered)}")
    print(f"arms with events  : {n_events} of {len(arms)}")
    print(f"panel genes        : {len(genes)}")
    print(f"ploidy            : {params['ploidy']}")
    print(f"min |log2|        : {threshold:.3f}"
          + (f" (noise sigma {noise_sigma:.4f})" if noise_sigma is not None
             else " (floor; too few baseline segments to estimate noise)"))
    print(f"below threshold   : {n_below} segments")
    print(f"baseline CN       : {baseline}"
          + ("" if baseline_fraction != baseline_fraction
             else f" ({baseline_fraction:.0%} of autosomal span)"))
    print(f"purity (flow)     : {args.purity if args.purity else 'not supplied'}")
    print(f"chromosome count  : {ploidy_calls['estimated_chromosome_count']}")
    print(f"ploidy class      : {ploidy_calls['ploidy_class']}"
          + (f" ({ploidy_calls['n_canonical_trisomies']} canonical trisomies)"
             if ploidy_calls["canonical_trisomies"] else ""))
    if ploidy_calls["trisomies"]:
        print(f"at three copies   : {', '.join(ploidy_calls['trisomies'])}")
    print(f"altered cytobands : {len(cytobands)}")
    limit = detection_limit(sigma_bin, bin_size, 5_000_000, baseline)
    print(f"detection limit   : "
          + (f"{limit:.1%} of cells for a 5 Mb change" if limit is not None
             else "not available"))
    print(f"reads per bin     : "
          + (f"{median_reads} median" if median_reads is not None
             else "no WIG supplied"))
    if median_reads is not None:
        print(f"overdispersion    : "
              f"{overdispersion_corrected}x corrected, {overdispersion}x raw "
              f"counts")
    for warning in warnings:
        print(f"WARNING           : {warning}")
    print(f"written           : {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
