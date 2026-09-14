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
    print(f"written           : {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
