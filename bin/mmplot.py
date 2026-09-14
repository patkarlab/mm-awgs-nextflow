#!/usr/bin/env python3
"""
mmplot.py

Copy number and allele fraction plots from mmkaryo.py and mmbaf.py output,
laid out to match the house style already used elsewhere in the lab: a depth
panel, a BAF panel with points coloured by call class, and a copy number panel
carrying total and minor allele copy number as separate traces, over a targets
track and a banded cytogram.

Modes:

  --genome   whole genome in one row, three panels plus a cytogram strip
  (default)  per-chromosome grid, each cell a compact three-panel stack
  --chr N    one chromosome full width, with gene highlight bands, a panel
             targets track and a labelled cytogram

Call classes in the BAF panel follow the segment each point falls in:

  DEL           total copy number below 2
  CNLOH         two copies with a minor allele of zero
  GAIN          more than two copies
  imbalance     allelically imbalanced without meeting the above
  balanced het  minor allele at half of total
  homozygous    no usable heterozygous signal

Usage:
  mmplot.py <prefix>.bins.tsv <prefix>.segments.baf.tsv <out.png> --purity P
            [--cytobands BED] [--genes BED] [--targets BED] [--chr chr16]

Dependencies: numpy, matplotlib.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from collections import OrderedDict

import numpy

CANONICAL = ["chr%d" % (i + 1) for i in range(22)] + ["chrX", "chrY"]

T2T_CENTROMERES = {
    "chr1": (121405145, 144008994), "chr2": (87403406, 98345164),
    "chr3": (90339372, 96498140), "chr4": (49050159, 57435402),
    "chr5": (46161727, 51082874), "chr6": (57216843, 62665581),
    "chr7": (55900922, 66709007), "chr8": (43506737, 48466602),
    "chr9": (40119875, 81269165), "chr10": (38527422, 43443011),
    "chr11": (48819609, 55216869), "chr12": (30997656, 38541929),
    "chr16": (32625033, 52263389), "chr17": (20527834, 28139193),
    "chr18": (11350790, 21125454), "chr19": (19876180, 30396484),
    "chr20": (25835813, 37114981),
}

GIEMSA = {"gneg": "#FFFFFF", "gpos25": "#D0D0D0", "gpos33": "#C0C0C0",
          "gpos50": "#9A9A9A", "gpos66": "#7A7A7A", "gpos75": "#5E5E5E",
          "gpos100": "#1C1C1C", "gpos": "#1C1C1C", "acen": "#B33A3A",
          "gvar": "#93ACCB", "stalk": "#6E93B8"}

CALL_COLOUR = OrderedDict([
    ("DEL", "#D1495B"),
    ("CNLOH", "#2E76BC"),
    ("GAIN", "#2E9E5B"),
    ("imbalance", "#8A5FBF"),
    ("balanced het", "#8A9098"),
    ("homozygous", "#DCDFE3"),
])

CN_TOTAL_COLOUR = "#1F4E79"
CN_MINOR_COLOUR = "#E07B22"

# Distinct fills for gene highlight bands, cycled.
HIGHLIGHT = ["#CFE2F3", "#D9EAD3", "#FCE5CD", "#EAD1DC", "#FFF2CC",
             "#D0E0E3", "#E6D5F2", "#F4CCCC"]


# ---------------------------------------------------------------------------
# Conversions
# ---------------------------------------------------------------------------

def total_copy_number(log2_value, purity, ploidy, reference_copies=2):
    if not math.isfinite(log2_value):
        return float("nan")
    denominator = purity * ploidy + (1.0 - purity) * 2.0
    return ((2.0 ** log2_value) * denominator
            - reference_copies * (1.0 - purity)) / purity


def minor_copy_number(theta, total, purity):
    if not (math.isfinite(theta) and math.isfinite(total)):
        return float("nan")
    denominator = purity * total + 2.0 * (1.0 - purity)
    return (theta * denominator - (1.0 - purity)) / purity


def call_class(total, minor, theta):
    """Classify a segment the way the BAF legend reads."""
    if not math.isfinite(total):
        return "balanced het"
    t = int(round(total))
    if math.isfinite(theta) and theta < 0.05:
        if t <= 0:
            return "homozygous"
    if t < 2:
        return "DEL"
    m = int(round(max(minor, 0.0))) if math.isfinite(minor) else None
    if t == 2:
        if m == 0:
            return "CNLOH"
        return "balanced het"
    if m is not None and m * 2 != t:
        return "GAIN"
    return "GAIN" if t > 2 else "imbalance"


# ---------------------------------------------------------------------------
# Readers
# ---------------------------------------------------------------------------

def _read_table(path, needed):
    if not path or not os.path.exists(path):
        return None
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: i for i, name in enumerate(header)}
        for name in needed:
            if name not in index:
                return None
        for line in handle:
            rows.append(line.rstrip("\n").split("\t"))
    return index, rows


def read_bins(path):
    parsed = _read_table(path, ["chromosome", "start", "log2_ratio"])
    if parsed is None:
        return {}
    index, rows = parsed
    by_chrom = {}
    for fields in rows:
        if fields[index["log2_ratio"]] == "NA":
            continue
        entry = by_chrom.setdefault(fields[index["chromosome"]], ([], []))
        entry[0].append(int(fields[index["start"]]))
        entry[1].append(float(fields[index["log2_ratio"]]))
    return {c: (numpy.array(v[0], dtype=numpy.int64), numpy.array(v[1]))
            for c, v in by_chrom.items()}


def read_windows(path):
    parsed = _read_table(path, ["chromosome", "start", "end", "theta"])
    if parsed is None:
        return {}
    index, rows = parsed
    by_chrom = {}
    for fields in rows:
        if fields[index["theta"]] == "NA":
            continue
        entry = by_chrom.setdefault(fields[index["chromosome"]], ([], []))
        entry[0].append((int(fields[index["start"]])
                         + int(fields[index["end"]])) // 2)
        entry[1].append(float(fields[index["theta"]]))
    return {c: (numpy.array(v[0], dtype=numpy.int64), numpy.array(v[1]))
            for c, v in by_chrom.items()}


def read_segments(path):
    parsed = _read_table(path, ["chromosome", "start", "end", "log2_ratio"])
    if parsed is None:
        raise SystemExit("ERROR: cannot read segments from %s" % path)
    index, rows = parsed
    segments = []
    for fields in rows:

        def get(name):
            if name not in index:
                return float("nan")
            raw = fields[index[name]]
            if raw in ("", "NA"):
                return float("nan")
            try:
                return float(raw)
            except ValueError:
                return float("nan")

        segments.append({"chrom": fields[index["chromosome"]],
                         "start": int(fields[index["start"]]),
                         "end": int(fields[index["end"]]),
                         "log2": get("log2_ratio"),
                         "theta": get("theta")})
    return segments


def read_bed(path, name_column=3):
    """Generic BED reader returning {chrom: [(start, end, name, extra), ...]}."""
    if not path or not os.path.exists(path):
        return {}
    by_chrom = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            name = fields[name_column] if len(fields) > name_column else ""
            extra = fields[4] if len(fields) > 4 else ""
            by_chrom.setdefault(chrom, []).append((start, end, name, extra))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


# ---------------------------------------------------------------------------
# Drawing helpers
# ---------------------------------------------------------------------------

def segment_lookup(positions, segments_here, key):
    """Value of `key` from the segment containing each position."""
    result = numpy.full(positions.shape[0], numpy.nan)
    objects = numpy.full(positions.shape[0], None, dtype=object)
    if not segments_here:
        return result, objects
    starts = numpy.array([s["start"] for s in segments_here])
    order = numpy.argsort(starts)
    starts = starts[order]
    ordered = [segments_here[i] for i in order]
    slot = numpy.searchsorted(starts, positions, side="right") - 1
    for i in numpy.nonzero((slot >= 0) & (slot < len(ordered)))[0]:
        segment = ordered[int(slot[i])]
        if positions[i] < segment["end"]:
            objects[i] = segment
            value = segment.get(key, numpy.nan)
            result[i] = value if isinstance(value, float) else numpy.nan
    return result, objects


def step_trace(axis, segments_here, key, colour, scale, linewidth=2.2,
               offset=0.0):
    """Segment values as a step trace, the way PURPLE draws copy number."""
    for s in sorted(segments_here, key=lambda x: x["start"]):
        value = s.get(key, float("nan"))
        if not math.isfinite(value):
            continue
        axis.plot([offset + s["start"] / scale, offset + s["end"] / scale],
                  [value] * 2, color=colour, linewidth=linewidth,
                  solid_capstyle="butt", zorder=4)


def draw_cytogram(axis, chrom, length, cytobands, label=False, offset=0.0,
                  scale=1e6):
    from matplotlib import patches
    bottom, height = 0.25, 0.50
    bands = cytobands.get(chrom, [])
    if bands:
        for start, end, name, stain in bands:
            axis.add_patch(patches.Rectangle(
                (offset + start / scale, bottom), (end - start) / scale, height,
                facecolor=GIEMSA.get(stain, "#FFFFFF"), edgecolor="0.45",
                linewidth=0.25, zorder=2))
            if label and (end - start) / scale > length / scale * 0.030:
                dark = stain in ("gpos100", "gpos75", "gpos", "acen")
                axis.text(offset + (start + end) / 2 / scale,
                          bottom + height / 2, name, ha="center", va="center",
                          fontsize=5.5, rotation=90,
                          color="white" if dark else "0.15", zorder=4)
    else:
        if chrom in T2T_CENTROMERES:
            a, b = T2T_CENTROMERES[chrom]
            arms = [(0, a), (b, length)]
        else:
            arms = [(0, length)]
        for a, b in arms:
            axis.add_patch(patches.Rectangle(
                (offset + a / scale, bottom), (b - a) / scale, height,
                facecolor="#E4E7EA", edgecolor="0.45", linewidth=0.6, zorder=2))
        if chrom in T2T_CENTROMERES:
            a, b = T2T_CENTROMERES[chrom]
            axis.add_patch(patches.Rectangle(
                (offset + a / scale, bottom + height * 0.2),
                (b - a) / scale, height * 0.6, facecolor="#B33A3A",
                edgecolor="none", zorder=3))
    axis.add_patch(patches.Rectangle(
        (offset, bottom), length / scale, height, facecolor="none",
        edgecolor="0.35", linewidth=0.7, zorder=5))


def style_track(axis):
    for spine in axis.spines.values():
        spine.set_visible(False)
    axis.set_yticks([])
    axis.set_ylim(0, 1)


def draw_panels(axes, chrom, length, bins, windows, segments, purity, ploidy,
                max_cn, point_size, genes=None, targets=None, cytobands=None,
                label_bands=False, compact=False):
    """Fill the depth, BAF, copy number, targets and cytogram axes."""
    ax_depth, ax_baf, ax_cn, ax_targets, ax_ideo = axes
    reference = 1 if chrom in ("chrX", "chrY") else 2
    here = [s for s in segments if s["chrom"] == chrom]
    scale = 1e6
    limit = length / scale

    # gene highlight bands behind everything
    if genes and not compact:
        for i, (start, end, name, _extra) in enumerate(genes.get(chrom, [])):
            colour = HIGHLIGHT[i % len(HIGHLIGHT)]
            for axis in (ax_depth, ax_baf, ax_cn):
                if axis is not None:
                    axis.axvspan(start / scale, end / scale, color=colour,
                                 alpha=0.55, zorder=0, linewidth=0)
            if ax_depth is not None and name:
                ax_depth.annotate(name, ((start + end) / 2 / scale, 1.02),
                                  xycoords=("data", "axes fraction"),
                                  ha="center", va="bottom", fontsize=8,
                                  fontweight="bold", color="0.2")

    for axis in [a for a in axes if a is not None]:
        if chrom in T2T_CENTROMERES and not (genes and not compact):
            a, b = T2T_CENTROMERES[chrom]
            axis.axvspan(a / scale, b / scale, color="0.95", zorder=0,
                         linewidth=0)
        axis.set_xlim(0, limit)

    # --- depth ---
    if chrom in bins:
        positions, log2 = bins[chrom]
        _v, objects = segment_lookup(positions, here, "log2")
        colours = [CALL_COLOUR[o["call"]] if o is not None
                   else CALL_COLOUR["balanced het"] for o in objects]
        ax_depth.scatter(positions / scale, log2, s=point_size, c=colours,
                         alpha=0.40, linewidths=0, zorder=2)
    step_trace(ax_depth, here, "log2", "black", scale, linewidth=2.4)
    ax_depth.axhline(0.0, color="0.45", linewidth=0.8, linestyle="--", zorder=1)
    ax_depth.set_ylim(-1.5, 1.5)
    ax_depth.tick_params(labelbottom=False, labelsize=8)

    # --- BAF ---
    if chrom in windows and chrom not in ("chrX", "chrY"):
        positions, theta = windows[chrom]
        _v, objects = segment_lookup(positions, here, "theta")
        colours = [CALL_COLOUR[o["call"]] if o is not None
                   else CALL_COLOUR["balanced het"] for o in objects]
        for values in (theta, 1.0 - theta):
            ax_baf.scatter(positions / scale, values, s=point_size * 1.5,
                           c=colours, alpha=0.60, linewidths=0, zorder=2)
    ax_baf.axhline(0.5, color="0.45", linewidth=0.8, zorder=1)
    ax_baf.set_ylim(0.0, 1.0)
    ax_baf.set_yticks([0.0, 0.5, 1.0])
    ax_baf.tick_params(labelbottom=False, labelsize=8)

    # --- copy number: total and minor allele ---
    step_trace(ax_cn, here, "total", CN_TOTAL_COLOUR, scale, linewidth=2.6)
    step_trace(ax_cn, here, "minor", CN_MINOR_COLOUR, scale, linewidth=2.2)
    for level in range(0, int(max_cn) + 1):
        ax_cn.axhline(level, color="0.93", linewidth=0.6, zorder=1)
    ax_cn.axhline(reference, color="0.55", linewidth=0.8, linestyle="--",
                  zorder=1)
    ax_cn.set_ylim(-0.4, max_cn)
    ax_cn.set_yticks(range(0, int(max_cn) + 1))
    ax_cn.tick_params(labelbottom=(ax_targets is None and ax_ideo is None),
                      labelsize=8)

    # --- panel targets ---
    if ax_targets is not None:
        from matplotlib import patches
        style_track(ax_targets)
        for start, end, name, _extra in (targets or {}).get(chrom, []):
            ax_targets.add_patch(patches.Rectangle(
                (start / scale, 0.30), max((end - start) / scale, limit * 0.0015),
                0.42, facecolor="#1F3B73", edgecolor="none", zorder=2))
        ax_targets.tick_params(labelbottom=False)

    # --- cytogram ---
    if ax_ideo is not None:
        style_track(ax_ideo)
        draw_cytogram(ax_ideo, chrom, length, cytobands or {}, label_bands)
        ax_ideo.tick_params(labelbottom=True, labelsize=8)


def main():
    parser = argparse.ArgumentParser(
        description="Copy number and allele fraction plots from "
                    "mmkaryo/mmbaf output")
    parser.add_argument("bins")
    parser.add_argument("segments")
    parser.add_argument("output")
    parser.add_argument("--purity", type=float, required=True,
                        help="Tumour purity, 0-1. Required: copy number cannot "
                             "be computed without it.")
    parser.add_argument("--ploidy", type=float, default=2.0)
    parser.add_argument("--chr", dest="chrom", default=None)
    parser.add_argument("--genome", action="store_true")
    parser.add_argument("--windows", default=None)
    parser.add_argument("--cytobands", default=None,
                        help="Cytoband BED (chrom, start, end, band, stain)")
    parser.add_argument("--genes", default=None,
                        help="Gene model BED, drawn as highlight bands with "
                             "labels in single-chromosome mode")
    parser.add_argument("--targets", default=None,
                        help="Panel BED, drawn as a targets track. Worth "
                             "showing even though these windows are masked "
                             "from the analysis, because the track marks where "
                             "the copy number calls are blind.")
    parser.add_argument("--max-cn", type=float, default=6.0)
    parser.add_argument("--columns", type=int, default=4)
    parser.add_argument("--title", default=None)
    parser.add_argument("--subtitle", default=None,
                        help="Extra text for the title line, e.g. a summary of "
                             "arm-level calls")
    args = parser.parse_args()

    try:
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib import pyplot, gridspec, lines
    except ImportError:
        raise SystemExit("ERROR: matplotlib is required.")

    rho, psi = args.purity, args.ploidy
    windows_path = args.windows or args.segments.replace(
        ".segments.baf.tsv", ".baf_windows.tsv")

    segments = read_segments(args.segments)
    bins = read_bins(args.bins)
    windows = read_windows(windows_path)
    cytobands = read_bed(args.cytobands)
    genes = read_bed(args.genes)
    targets = read_bed(args.targets)

    for s in segments:
        reference = 1 if s["chrom"] in ("chrX", "chrY") else 2
        s["total"] = total_copy_number(s["log2"], rho, psi, reference)
        s["minor"] = minor_copy_number(s["theta"], s["total"], rho)
        s["call"] = call_class(s["total"], s["minor"], s["theta"])

    lengths = {}
    for s in segments:
        lengths[s["chrom"]] = max(lengths.get(s["chrom"], 0), s["end"])
    for chrom, (positions, _v) in bins.items():
        lengths[chrom] = max(lengths.get(chrom, 0), int(positions.max()))

    title = args.title or os.path.basename(args.segments).split(".")[0]
    header = "%s    purity %.2f, ploidy %.2f" % (title, rho, psi)
    if args.subtitle:
        header += "   |   " + args.subtitle

    baf_handles = [lines.Line2D([], [], marker="o", linestyle="none",
                                markersize=5, color=c, label=k)
                   for k, c in CALL_COLOUR.items()]
    cn_handles = [
        lines.Line2D([], [], color=CN_TOTAL_COLOUR, lw=3, label="total CN"),
        lines.Line2D([], [], color=CN_MINOR_COLOUR, lw=3, label="minor allele CN"),
    ]

    # ---- single chromosome ------------------------------------------------
    if args.chrom:
        chrom = args.chrom
        if chrom not in lengths:
            raise SystemExit("ERROR: no data for %s" % chrom)
        figure = pyplot.figure(figsize=(17, 10))
        grid = gridspec.GridSpec(5, 1, height_ratios=[3, 2.2, 2.4, 0.45, 0.85],
                                 hspace=0.10)
        axes = [figure.add_subplot(grid[0])]
        for i in range(1, 5):
            axes.append(figure.add_subplot(grid[i], sharex=axes[0]))
        draw_panels(axes, chrom, lengths[chrom], bins, windows, segments,
                    rho, psi, args.max_cn, point_size=13.0, genes=genes,
                    targets=targets, cytobands=cytobands, label_bands=True)
        axes[0].set_ylabel("log2 depth ratio")
        axes[1].set_ylabel("BAF\nallele fraction")
        axes[2].set_ylabel("copy number")
        axes[3].set_ylabel("targets", fontsize=8, rotation=0, ha="right",
                           va="center")
        axes[4].set_ylabel("cytoband", fontsize=8, rotation=0, ha="right",
                           va="center")
        axes[4].set_xlabel("%s position (Mb)" % chrom)
        axes[0].set_title("%s  --  %s" % (chrom, header), fontsize=13, pad=22)
        axes[1].legend(handles=baf_handles, fontsize=7, ncol=6,
                       loc="upper right", framealpha=0.92)
        axes[2].legend(handles=cn_handles, fontsize=8, ncol=2,
                       loc="upper right", framealpha=0.92)
        figure.savefig(args.output, bbox_inches="tight", dpi=130)
        sys.stderr.write("wrote %s\n" % args.output)
        return

    order = [c for c in CANONICAL if c in lengths]

    # ---- genome wide ------------------------------------------------------
    if args.genome:
        offsets, running = OrderedDict(), 0
        for name in order:
            offsets[name] = running
            running += lengths[name]
        scale = 1e6
        figure = pyplot.figure(figsize=(21, 10))
        grid = gridspec.GridSpec(4, 1, height_ratios=[3, 2.2, 2.4, 0.7],
                                 hspace=0.10)
        ax_depth = figure.add_subplot(grid[0])
        ax_baf = figure.add_subplot(grid[1], sharex=ax_depth)
        ax_cn = figure.add_subplot(grid[2], sharex=ax_depth)
        ax_ideo = figure.add_subplot(grid[3], sharex=ax_depth)

        for name in order:
            shift = offsets[name] / scale
            here = [s for s in segments if s["chrom"] == name]
            if name in bins:
                positions, log2 = bins[name]
                _v, objects = segment_lookup(positions, here, "log2")
                colours = [CALL_COLOUR[o["call"]] if o is not None
                           else CALL_COLOUR["balanced het"] for o in objects]
                ax_depth.scatter(shift + positions / scale, log2, s=1.5,
                                 c=colours, alpha=0.35, linewidths=0, zorder=2)
            if name in windows and name not in ("chrX", "chrY"):
                positions, theta = windows[name]
                _v, objects = segment_lookup(positions, here, "theta")
                colours = [CALL_COLOUR[o["call"]] if o is not None
                           else CALL_COLOUR["balanced het"] for o in objects]
                for values in (theta, 1.0 - theta):
                    ax_baf.scatter(shift + positions / scale, values, s=3.0,
                                   c=colours, alpha=0.5, linewidths=0, zorder=2)
            step_trace(ax_depth, here, "log2", "black", scale, 2.2, shift)
            step_trace(ax_cn, here, "total", CN_TOTAL_COLOUR, scale, 2.4, shift)
            step_trace(ax_cn, here, "minor", CN_MINOR_COLOUR, scale, 2.0, shift)
            draw_cytogram(ax_ideo, name, lengths[name], cytobands,
                          label=False, offset=shift)

        ax_depth.axhline(0.0, color="0.45", linewidth=0.8, linestyle="--")
        ax_depth.set_ylim(-1.5, 1.5)
        ax_depth.set_ylabel("depth\nlog2 copy ratio")
        ax_baf.axhline(0.5, color="0.45", linewidth=0.8)
        ax_baf.set_ylim(0, 1)
        ax_baf.set_yticks([0.0, 0.5, 1.0])
        ax_baf.set_ylabel("BAF\nallele fraction")
        for level in range(0, int(args.max_cn) + 1):
            ax_cn.axhline(level, color="0.93", linewidth=0.6, zorder=1)
        ax_cn.axhline(2, color="0.55", linewidth=0.8, linestyle="--", zorder=1)
        ax_cn.set_ylim(-0.4, args.max_cn)
        ax_cn.set_yticks(range(0, int(args.max_cn) + 1))
        ax_cn.set_ylabel("copy number")
        style_track(ax_ideo)

        for axis in (ax_depth, ax_baf, ax_cn, ax_ideo):
            for name in order[:-1]:
                axis.axvline((offsets[name] + lengths[name]) / scale,
                             color="0.88", linewidth=0.6, zorder=1)
            axis.set_xlim(0, running / scale)
            axis.tick_params(labelbottom=False, labelsize=8)
        ax_ideo.tick_params(labelbottom=True, labelsize=9)
        ax_ideo.set_xticks([(offsets[n] + lengths[n] / 2) / scale
                            for n in order])
        ax_ideo.set_xticklabels([n[3:] for n in order], fontsize=9)
        ax_ideo.set_xlabel("chromosome")
        ax_depth.set_title(header, fontsize=13)
        ax_baf.legend(handles=baf_handles, fontsize=7, ncol=6,
                      loc="upper right", framealpha=0.92)
        ax_cn.legend(handles=cn_handles, fontsize=8, ncol=2,
                     loc="upper right", framealpha=0.92)
        figure.savefig(args.output, bbox_inches="tight", dpi=120)
        sys.stderr.write("wrote %s\n" % args.output)
        return

    # ---- grid -------------------------------------------------------------
    columns = max(1, args.columns)
    rows = int(math.ceil(len(order) / float(columns)))
    figure = pyplot.figure(figsize=(5.4 * columns, 4.4 * rows))
    outer = gridspec.GridSpec(rows, columns, hspace=0.48, wspace=0.22)

    for position, chrom in enumerate(order):
        cell = gridspec.GridSpecFromSubplotSpec(
            4, 1, subplot_spec=outer[position],
            height_ratios=[2.4, 1.8, 2.0, 0.5], hspace=0.10)
        ax_depth = figure.add_subplot(cell[0])
        ax_baf = figure.add_subplot(cell[1], sharex=ax_depth)
        ax_cn = figure.add_subplot(cell[2], sharex=ax_depth)
        ax_ideo = figure.add_subplot(cell[3], sharex=ax_depth)
        draw_panels([ax_depth, ax_baf, ax_cn, None, ax_ideo], chrom,
                    lengths[chrom], bins, windows, segments, rho, psi,
                    args.max_cn, point_size=2.4, cytobands=cytobands,
                    compact=True)
        ax_depth.set_title(chrom, fontsize=11, pad=4)
        if position % columns == 0:
            ax_depth.set_ylabel("log2", fontsize=9)
            ax_baf.set_ylabel("BAF", fontsize=9)
            ax_cn.set_ylabel("CN", fontsize=9)

    figure.suptitle(header, fontsize=15, y=0.998)
    figure.legend(handles=baf_handles + cn_handles, fontsize=10, ncol=8,
                  loc="lower center", bbox_to_anchor=(0.5, -0.010),
                  framealpha=0.92)
    figure.savefig(args.output, bbox_inches="tight", dpi=110)
    sys.stderr.write("wrote %s\n" % args.output)


if __name__ == "__main__":
    main()
