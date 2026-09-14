#!/usr/bin/env python3
"""
compare_sv_callers.py

Compare structural variant calls across Sniffles, CuteSV, Severus and SAVANA
for one sample, clustering junctions that the callers describe at slightly
different coordinates and reporting which callers found each one.

Written to answer a specific question before SAVANA is wired into the SURVIVOR
ensemble: what does adding roughly 2,500 unfiltered calls per sample do to the
merged output, and are the additional calls junctions the other three missed or
artifacts of on-target alignment noise.

No expected translocation, FISH finding or known breakpoint is encoded
anywhere. The script reports every junction it finds and annotates each against
the panel BED and the Ig segment BED; which of them matter is for the reader to
decide.

Junction normalisation
----------------------
Callers describe the same junction differently. Sniffles and CuteSV emit TRA
records with CHR2 and END in INFO; SAVANA and Severus emit paired BND records
whose ALT field carries the mate coordinate. Both forms reduce to an unordered
pair of breakends, which is what gets clustered.

Two junctions are treated as the same event when both breakends fall within
--tolerance of each other on matching chromosomes. Default 500 bp, which is
generous for long reads but accommodates the coordinate differences seen
between callers on the same junction.

Artifact flagging
-----------------
Intrachromosomal junctions shorter than --min-intra (default 1000 bp) are
flagged rather than dropped. These dominate SAVANA's unfiltered output: on one
sample every call in the annotated set was a self-junction inside a single
panel gene body, 33 to 369 bp across, which is local alignment noise at high
on-target depth rather than structural variation. Flagging rather than
discarding keeps the count visible.

Support fields, per caller:
    Sniffles   INFO/SUPPORT
    CuteSV     INFO/RE
    Severus    FORMAT/DV
    SAVANA     INFO/TUMOUR_READ_SUPPORT

SAVANA's CLASS field is carried through as an annotation column. It is not used
to filter: on this data the somatic model classifies low-support translocations
as noise, including ones independently confirmed by another caller.

Dependencies: standard library only.
"""

import argparse
import csv
import gzip
import os
import re
import sys
from collections import defaultdict

BND_ALT = re.compile(r'[\[\]]([^\[\]:]+):(\d+)[\[\]]')

SUPPORT_FIELDS = {
    "sniffles": ("INFO", "SUPPORT"),
    "cutesv": ("INFO", "RE"),
    "severus": ("FORMAT", "DV"),
    "savana": ("INFO", "TUMOUR_READ_SUPPORT"),
}


def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def parse_info(field):
    info = {}
    for item in field.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        elif item:
            info[item] = True
    return info


def get_support(caller, info, fields):
    """Pull the caller's read-support value, or None if absent."""
    source, key = SUPPORT_FIELDS.get(caller, ("INFO", None))
    if key is None:
        return None

    if source == "INFO":
        value = info.get(key)
    else:
        if len(fields) < 10:
            return None
        format_keys = fields[8].split(":")
        sample_values = fields[9].split(":")
        if key not in format_keys:
            return None
        index = format_keys.index(key)
        value = sample_values[index] if index < len(sample_values) else None

    if value in (None, True, "", "."):
        return None
    try:
        return int(str(value).split(",")[0])
    except ValueError:
        return None


def mate_from_alt(alt):
    match = BND_ALT.search(alt)
    if not match:
        return None, None
    return match.group(1), int(match.group(2))


def read_vcf(path, caller):
    """Return a list of junction dicts from one caller's VCF."""
    records = []
    seen_mates = set()

    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom_a, pos_a = fields[0], int(fields[1])
            record_id, alt, filt = fields[2], fields[4], fields[6]
            info = parse_info(fields[7])
            svtype = info.get("SVTYPE", "")

            chrom_b = pos_b = None

            if svtype == "BND" or "[" in alt or "]" in alt:
                chrom_b, pos_b = mate_from_alt(alt)
                # Paired BND records describe one junction twice. Keep the
                # first and skip its mate.
                mate = info.get("MATEID")
                if mate and mate in seen_mates:
                    continue
                if mate:
                    seen_mates.add(record_id)
            elif svtype in ("TRA", "BND_TRA"):
                chrom_b = info.get("CHR2")
                end = info.get("END")
                pos_b = int(end) if end and str(end).lstrip("-").isdigit() else None
            else:
                # Intrachromosomal typed event (DEL/DUP/INV/INS).
                chrom_b = info.get("CHR2", chrom_a)
                end = info.get("END")
                pos_b = int(end) if end and str(end).lstrip("-").isdigit() else None

            if chrom_b is None or pos_b is None:
                continue

            records.append({
                "caller": caller,
                "id": record_id,
                "svtype": svtype or "NA",
                "filter": filt,
                "chrom_a": chrom_a, "pos_a": pos_a,
                "chrom_b": chrom_b, "pos_b": pos_b,
                "support": get_support(caller, info, fields),
                "savana_class": info.get("CLASS", ""),
            })

    return records


def canonical(record):
    """Order the two breakends so the same junction always keys identically."""
    a = (record["chrom_a"], record["pos_a"])
    b = (record["chrom_b"], record["pos_b"])
    return (a, b) if a <= b else (b, a)


def cluster(records, tolerance):
    """Group records describing the same junction. Simple greedy clustering."""
    clusters = []
    for record in sorted(records, key=canonical):
        (chrom_a, pos_a), (chrom_b, pos_b) = canonical(record)
        placed = False
        for group in clusters:
            if group["chrom_a"] != chrom_a or group["chrom_b"] != chrom_b:
                continue
            if (abs(group["pos_a"] - pos_a) <= tolerance and
                    abs(group["pos_b"] - pos_b) <= tolerance):
                group["records"].append(record)
                placed = True
                break
        if not placed:
            clusters.append({
                "chrom_a": chrom_a, "pos_a": pos_a,
                "chrom_b": chrom_b, "pos_b": pos_b,
                "records": [record],
            })
    return clusters


def load_bed(path):
    intervals = defaultdict(list)
    if not path or not os.path.isfile(path):
        return dict(intervals)
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split()
            name = parts[3] if len(parts) > 3 else ""
            intervals[parts[0]].append((int(parts[1]), int(parts[2]), name))
    for chrom in intervals:
        intervals[chrom].sort()
    return dict(intervals)


def annotate_position(intervals, chrom, pos):
    """Return (name_if_inside, distance_to_nearest_boundary)."""
    if chrom not in intervals:
        return "", ""
    best_name, best_distance = "", None
    for start, end, name in intervals[chrom]:
        if start <= pos < end:
            best_name = name or "in_panel"
            distance = min(pos - start, end - pos)
        else:
            distance = start - pos if pos < start else pos - end
        if best_distance is None or distance < best_distance:
            best_distance = distance
    return best_name, ("" if best_distance is None else best_distance)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sniffles")
    parser.add_argument("--cutesv")
    parser.add_argument("--severus")
    parser.add_argument("--savana")
    parser.add_argument("--panel-bed", help="capture BED for this sample")
    parser.add_argument("--ig-bed", help="Ig segment BED (assets/ig_segments_t2t.bed)")
    parser.add_argument("--out", required=True, help="output TSV")
    parser.add_argument("--tolerance", type=int, default=500,
                        help="breakend distance for two calls to be one junction")
    parser.add_argument("--min-intra", type=int, default=1000,
                        help="intrachromosomal junctions shorter than this are "
                             "flagged as likely local alignment artifacts")
    parser.add_argument("--interchromosomal-only", action="store_true",
                        help="restrict output to junctions between chromosomes")
    args = parser.parse_args()

    sources = [("sniffles", args.sniffles), ("cutesv", args.cutesv),
               ("severus", args.severus), ("savana", args.savana)]

    records = []
    present = []
    for caller, path in sources:
        if not path:
            continue
        if not os.path.isfile(path):
            print(f"WARNING: {caller} VCF not found: {path}", file=sys.stderr)
            continue
        found = read_vcf(path, caller)
        records.extend(found)
        present.append(caller)
        print(f"{caller:9} {len(found):7d} junctions from {os.path.basename(path)}")

    if not records:
        sys.exit("ERROR: no calls read from any caller")

    panel = load_bed(args.panel_bed)
    ig = load_bed(args.ig_bed)

    clusters = cluster(records, args.tolerance)

    rows = []
    for group in clusters:
        interchromosomal = group["chrom_a"] != group["chrom_b"]
        span = None if interchromosomal else abs(group["pos_b"] - group["pos_a"])

        if args.interchromosomal_only and not interchromosomal:
            continue

        callers = sorted({r["caller"] for r in group["records"]})
        support = {}
        for caller in present:
            values = [r["support"] for r in group["records"]
                      if r["caller"] == caller and r["support"] is not None]
            support[caller] = max(values) if values else ""

        savana_classes = sorted({r["savana_class"] for r in group["records"]
                                 if r["caller"] == "savana" and r["savana_class"]})

        panel_a, dist_a = annotate_position(panel, group["chrom_a"], group["pos_a"])
        panel_b, dist_b = annotate_position(panel, group["chrom_b"], group["pos_b"])
        ig_a, _ = annotate_position(ig, group["chrom_a"], group["pos_a"])
        ig_b, _ = annotate_position(ig, group["chrom_b"], group["pos_b"])

        flags = []
        if not interchromosomal and span is not None and span < args.min_intra:
            flags.append("short_intrachromosomal")
        if len(callers) == 1:
            flags.append(f"{callers[0]}_only")
        if ig_a or ig_b:
            flags.append("ig_breakend")

        rows.append({
            "sample": args.sample,
            "chrom_a": group["chrom_a"], "pos_a": group["pos_a"],
            "chrom_b": group["chrom_b"], "pos_b": group["pos_b"],
            "interchromosomal": "yes" if interchromosomal else "no",
            "span": "" if span is None else span,
            "n_callers": len(callers),
            "callers": ",".join(callers),
            **{f"support_{c}": support.get(c, "") for c in present},
            "savana_class": ",".join(savana_classes),
            "panel_region_a": panel_a, "panel_dist_a": dist_a,
            "panel_region_b": panel_b, "panel_dist_b": dist_b,
            "ig_region_a": ig_a, "ig_region_b": ig_b,
            "flags": ";".join(flags),
        })

    rows.sort(key=lambda r: (-r["n_callers"], r["chrom_a"], r["pos_a"]))

    fieldnames = (["sample", "chrom_a", "pos_a", "chrom_b", "pos_b",
                   "interchromosomal", "span", "n_callers", "callers"] +
                  [f"support_{c}" for c in present] +
                  ["savana_class", "panel_region_a", "panel_dist_a",
                   "panel_region_b", "panel_dist_b", "ig_region_a", "ig_region_b",
                   "flags"])

    with open(args.out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"\njunction clusters: {len(rows)}")

    by_count = defaultdict(int)
    for row in rows:
        by_count[row["n_callers"]] += 1
    for n in sorted(by_count, reverse=True):
        print(f"  found by {n} caller(s): {by_count[n]}")

    inter = [r for r in rows if r["interchromosomal"] == "yes"]
    print(f"\ninterchromosomal: {len(inter)}")
    ig_rows = [r for r in rows if r["ig_region_a"] or r["ig_region_b"]]
    print(f"with an Ig breakend: {len(ig_rows)}")
    short = [r for r in rows if "short_intrachromosomal" in r["flags"]]
    print(f"flagged short intrachromosomal: {len(short)}")

    for caller in present:
        only = [r for r in rows if r["callers"] == caller]
        print(f"{caller}-only: {len(only)}")

    print(f"\nwritten: {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
