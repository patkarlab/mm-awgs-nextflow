#!/usr/bin/env python3
"""
plot_baf_cn.py

Combined B-allele frequency and copy number visualisation for panel-based
adaptive sampling data.

Rationale
---------
Copy-neutral LOH is defined by allelic imbalance in the absence of copy number
change. BAF alone cannot distinguish it from hemizygous deletion, and copy
number alone cannot detect it at all. The two tracks must therefore be read
together on a shared genomic axis, which is what this script produces.

Three views are generated:

  1. Genome-wide, per sample. Copy number from ichorCNA bins on top, panel-window
     BAF below, sharing an x-axis. Panel windows are shaded so it is visually
     obvious where BAF data can and cannot exist. This is the orienting view and
     the one that guards against over-reading sparse off-target regions.

  2. Per panel region, per sample. Individual het sites plotted as points against
     position, so bimodality is directly visible rather than reduced to a summary
     statistic. Regions failing the assessability gate are drawn but explicitly
     annotated as unassessable.

  3. Cohort heatmap, samples against panel regions, coloured by BAF deflection.
     Deflection is shown rather than the categorical flag so that near-threshold
     regions remain visible. Unassessable cells are hatched rather than coloured,
     since any colour drawn from the scale would be read as a measured value.

Copy number is re-plotted from the ichorCNA .cna.seg output rather than reusing
the PDF ichorCNA emits. Sharing an x-axis exactly between the CN and BAF panels
is the entire point of the figure, and that alignment cannot be achieved by
annotating a fixed-layout PDF.

Dependencies
------------
  Python standard library plus matplotlib. No pandas.

Usage
-----
  python3 plot_baf_cn.py \
      --screen   cohort.baf_screen.tsv \
      --bed      aWGS_PCN_v7_t2t_chr.bed \
      --sample-map plot_inputs.tsv \
      --outdir   results/baf_cn

  plot_inputs.tsv is tab-separated with a header, columns:
      sample, vcf_path, cna_seg, params_txt
  cna_seg and params_txt may be given as NA when ichorCNA output is unavailable,
  in which case the BAF panel is drawn without a copy number track.
"""

import argparse
import glob
import gzip
import os
import statistics
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Patch


# ---------------------------------------------------------------------------
# Presentation constants
# ---------------------------------------------------------------------------

COLOUR_BALANCED = "#3B6FA0"
COLOUR_LOH = "#B03A2E"
COLOUR_EQUIVOCAL = "#D68910"
COLOUR_CN = "#4A4A4A"
COLOUR_PANEL_SHADE = "#DCE6F0"
COLOUR_UNASSESSABLE = "#BFBFBF"

FLAG_COLOURS = {
    "LOH_LIKELY": COLOUR_LOH,
    "EQUIVOCAL": COLOUR_EQUIVOCAL,
    "NO_LOH": COLOUR_BALANCED,
    "UNASSESSABLE": COLOUR_UNASSESSABLE,
}

# Chromosome ordering for genome-wide layout. Chromosomes absent from the data
# are skipped, so this list is safe to apply to either reference naming.
CHROM_ORDER = ["chr{0}".format(i) for i in range(1, 23)] + ["chrX", "chrY"]


# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------

def open_maybe_gzip(path):
    """Return a text-mode handle for a plain or gzip-compressed file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def read_bed(path):
    """Read a four-column panel BED into a list of region dicts."""
    regions = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            regions.append({
                "name": fields[3] if len(fields) > 3 else "{0}:{1}".format(fields[0], fields[1]),
                "chrom": fields[0],
                "start": int(fields[1]),
                "end": int(fields[2]),
            })
    return regions


def read_screen_table(path):
    """
    Read the output of baf_loh_screen.py.

    Returns a dict keyed by (sample, region) and the ordered list of sample
    identifiers encountered.
    """
    rows = {}
    samples = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name: position for position, name in enumerate(header)}
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue

            def value(name):
                raw = fields[index[name]] if name in index else "NA"
                return None if raw == "NA" else raw

            sample = value("sample")
            region = value("region")
            if sample is None or region is None:
                continue
            if sample not in samples:
                samples.append(sample)

            def as_float(name):
                raw = value(name)
                try:
                    return float(raw) if raw is not None else None
                except ValueError:
                    return None

            rows[(sample, region)] = {
                "chrom": value("chrom"),
                "start": int(value("start")) if value("start") else None,
                "end": int(value("end")) if value("end") else None,
                "n_het": int(value("n_het")) if value("n_het") else 0,
                "median_dp": as_float("median_dp"),
                "median_baf": as_float("median_baf"),
                "frac_central": as_float("frac_central"),
                "bimodality": as_float("bimodality"),
                "baf_deflection": as_float("baf_deflection"),
                "depletion_score": as_float("depletion_score"),
                "flag": value("flag") or "UNASSESSABLE",
            }
    return rows, samples


def read_ichor_params(path):
    """
    Extract tumour fraction and ploidy from an ichorCNA .params.txt file.

    The file is a small key/value table; only the two headline values are
    needed here, for figure annotation.
    """
    result = {"tumour_fraction": None, "ploidy": None}
    if not path or path == "NA" or not os.path.exists(path):
        return result
    with open(path) as handle:
        for line in handle:
            fields = [field.strip() for field in line.rstrip("\n").split("\t")]
            if len(fields) < 2:
                continue
            key = fields[0].lower()
            try:
                number = float(fields[1])
            except ValueError:
                continue
            if "tumor fraction" in key or "tumour fraction" in key:
                result["tumour_fraction"] = number
            elif "ploidy" in key:
                result["ploidy"] = number
    return result


def read_cna_seg(path):
    """
    Read bin-level copy number from an ichorCNA .cna.seg file.

    Returns a dict keyed by chromosome, each value a list of
    (position, log2_ratio) tuples. Column names vary slightly between ichorCNA
    versions, so the log ratio column is located by header inspection rather
    than by fixed index.
    """
    if not path or path == "NA" or not os.path.exists(path):
        return {}

    by_chrom = {}
    with open(path) as handle:
        header_line = handle.readline().rstrip("\n")
        header = [column.strip().lower() for column in header_line.split("\t")]

        def find_column(candidates):
            for candidate in candidates:
                for position, name in enumerate(header):
                    if name == candidate or name.endswith("." + candidate):
                        return position
            for candidate in candidates:
                for position, name in enumerate(header):
                    if candidate in name:
                        return position
            return None

        chrom_index = find_column(["chr", "chrom", "chromosome"])
        start_index = find_column(["start"])
        ratio_index = find_column(["logr", "log2_ratio_median", "median", "copy.number.log"])

        if chrom_index is None or start_index is None or ratio_index is None:
            sys.stderr.write(
                "Could not locate required columns in {0}; copy number track skipped\n".format(path)
            )
            return {}

        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(chrom_index, start_index, ratio_index):
                continue
            chrom = fields[chrom_index]
            if not chrom.startswith("chr"):
                chrom = "chr" + chrom
            try:
                start = int(float(fields[start_index]))
                ratio = float(fields[ratio_index])
            except ValueError:
                continue
            by_chrom.setdefault(chrom, []).append((start, ratio))

    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda item: item[0])
    return by_chrom


# ---------------------------------------------------------------------------
# Het site loading for point-level plots
#
# The screen script summarises sites; the per-region plots need the sites
# themselves. The parsing logic is intentionally kept simple and mirrors the
# screen's filters so the two views cannot disagree.
# ---------------------------------------------------------------------------

def collect_vcf_paths(target):
    """Resolve a VCF file or directory of per-chromosome VCFs to a path list."""
    if os.path.isfile(target):
        return [target]
    if os.path.isdir(target):
        paths = []
        for pattern in ("phased_*.vcf.gz", "phased_*.vcf"):
            paths.extend(glob.glob(os.path.join(target, pattern)))
        seen = set()
        unique = []
        for path in sorted(paths):
            stem = path[:-3] if path.endswith(".gz") else path
            if stem not in seen:
                seen.add(stem)
                unique.append(path)
        return unique
    return []


def load_het_sites(vcf_paths, min_site_depth):
    """
    Load het SNV sites as {chrom: [(pos, depth, minor_baf), ...]}.

    Filters match baf_loh_screen.py: PASS biallelic SNVs, diploid heterozygous
    genotype, depth at or above the minimum.
    """
    by_chrom = {}
    for vcf_path in vcf_paths:
        with open_maybe_gzip(vcf_path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 10:
                    continue
                ref, alt, filt = fields[3], fields[4], fields[6]
                if filt not in (".", "PASS", ""):
                    continue
                if "," in alt or len(ref) != 1 or len(alt) != 1:
                    continue

                keys = fields[8].split(":")
                values = fields[9].split(":")
                record = dict(zip(keys, values))

                genotype = record.get("GT", "./.")
                alleles = genotype.replace("|", "/").split("/")
                if len(alleles) != 2 or "." in alleles or alleles[0] == alleles[1]:
                    continue

                allele_depths = record.get("AD", "")
                depth = None
                minor = None
                if allele_depths and allele_depths != ".":
                    try:
                        counts = [int(x) for x in allele_depths.split(",") if x != "."]
                    except ValueError:
                        counts = []
                    if len(counts) >= 2 and (counts[0] + counts[1]) > 0:
                        depth = counts[0] + counts[1]
                        # Unfolded alternate-allele fraction. LOH presents as
                        # two lobes either side of 0.5; folding would collapse
                        # them together and hide the signature.
                        minor = counts[1] / float(depth)
                if depth is None:
                    continue
                if depth < min_site_depth:
                    continue

                try:
                    position = int(fields[1])
                except ValueError:
                    continue
                by_chrom.setdefault(fields[0], []).append((position, depth, minor))

    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda item: item[0])
    return by_chrom


# ---------------------------------------------------------------------------
# Genome-wide coordinate layout
# ---------------------------------------------------------------------------

def build_genome_layout(cna_by_chrom, regions):
    """
    Build cumulative genomic offsets so all chromosomes lay out on one axis.

    Chromosome extents are taken from the union of the ichorCNA bins and the
    panel regions, so the layout works whether or not ichorCNA output is
    present.
    """
    extents = {}
    for chrom, bins in cna_by_chrom.items():
        if bins:
            extents[chrom] = max(extents.get(chrom, 0), bins[-1][0])
    for region in regions:
        extents[region["chrom"]] = max(extents.get(region["chrom"], 0), region["end"])

    ordered = [chrom for chrom in CHROM_ORDER if chrom in extents]
    ordered.extend(sorted(chrom for chrom in extents if chrom not in CHROM_ORDER))

    offsets = {}
    cumulative = 0
    for chrom in ordered:
        offsets[chrom] = cumulative
        cumulative += extents[chrom]
    return offsets, ordered, cumulative


# ---------------------------------------------------------------------------
# Figure 1: genome-wide CN over BAF
# ---------------------------------------------------------------------------

def plot_genome_wide(sample, screen_rows, regions, het_by_chrom,
                     cna_by_chrom, params, outdir, min_sites):
    """
    Draw the genome-wide two-panel figure for one sample.

    Top panel is copy number, bottom panel is minor-allele BAF. Panel windows
    are shaded across both panels so that the reader can see at a glance which
    parts of the genome carry BAF information at all. Regions that failed the
    assessability gate are hatched rather than left blank, so a gap is never
    read as evidence of normality.
    """
    offsets, ordered_chroms, genome_length = build_genome_layout(cna_by_chrom, regions)
    if genome_length <= 0:
        sys.stderr.write("{0}: no coordinates available, skipping genome-wide plot\n".format(sample))
        return None

    figure, (cn_axis, baf_axis) = plt.subplots(
        2, 1, figsize=(16, 7), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.3], "hspace": 0.08},
    )

    # Shade panel windows on both panels.
    for region in regions:
        if region["chrom"] not in offsets:
            continue
        offset = offsets[region["chrom"]]
        row = screen_rows.get((sample, region["name"]))
        unassessable = (row is None or row["flag"] == "UNASSESSABLE")
        for axis in (cn_axis, baf_axis):
            axis.axvspan(
                offset + region["start"], offset + region["end"],
                color=COLOUR_PANEL_SHADE, zorder=0,
                hatch="///" if unassessable else None,
                edgecolor="#9AA5B1" if unassessable else None,
                linewidth=0.0,
            )

    # Copy number panel.
    if cna_by_chrom:
        x_values = []
        y_values = []
        for chrom in ordered_chroms:
            offset = offsets[chrom]
            for position, ratio in cna_by_chrom.get(chrom, []):
                x_values.append(offset + position)
                y_values.append(ratio)
        cn_axis.scatter(x_values, y_values, s=1.2, c=COLOUR_CN, alpha=0.35,
                        linewidths=0, rasterized=True, zorder=2)
        cn_axis.axhline(0.0, color="#666666", linewidth=0.9, linestyle="-", zorder=3)
        cn_axis.set_ylim(-2.0, 2.0)
    else:
        cn_axis.text(0.5, 0.5, "No ichorCNA copy number available",
                     transform=cn_axis.transAxes, ha="center", va="center",
                     fontsize=10, color="#888888")
        cn_axis.set_ylim(-2.0, 2.0)

    cn_axis.set_ylabel("Copy number\nlog2 ratio", fontsize=10)

    annotation = []
    if params.get("tumour_fraction") is not None:
        annotation.append("tumour fraction {0:.2f}".format(params["tumour_fraction"]))
    if params.get("ploidy") is not None:
        annotation.append("ploidy {0:.2f}".format(params["ploidy"]))
    title = "Sample {0} - copy number and minor-allele BAF".format(sample)
    if annotation:
        title += "   ({0})".format(", ".join(annotation))
    cn_axis.set_title(title, fontsize=12, pad=10)

    # BAF panel: individual het sites inside panel windows only.
    for region in regions:
        if region["chrom"] not in offsets:
            continue
        offset = offsets[region["chrom"]]
        row = screen_rows.get((sample, region["name"]))
        flag = row["flag"] if row else "UNASSESSABLE"
        if flag == "UNASSESSABLE":
            continue
        sites = [
            (position, minor)
            for position, _depth, minor in het_by_chrom.get(region["chrom"], [])
            if region["start"] <= position - 1 < region["end"]
        ]
        if not sites:
            continue
        baf_axis.scatter(
            [offset + position for position, _minor in sites],
            [minor for _position, minor in sites],
            s=2.0, c=FLAG_COLOURS.get(flag, COLOUR_BALANCED),
            alpha=0.45, linewidths=0, rasterized=True, zorder=2,
        )

    baf_axis.axhline(0.5, color="#666666", linewidth=0.9, zorder=3)
    baf_axis.axhspan(0.42, 0.58, color="#2E7D32", alpha=0.07, zorder=1)
    baf_axis.set_ylim(0.0, 1.0)
    baf_axis.set_ylabel("B-allele frequency", fontsize=10)

    # Chromosome boundaries and labels.
    for axis in (cn_axis, baf_axis):
        for chrom in ordered_chroms:
            axis.axvline(offsets[chrom], color="#CCCCCC", linewidth=0.6, zorder=1)
        axis.set_xlim(0, genome_length)
        axis.grid(axis="y", color="#EEEEEE", linewidth=0.6, zorder=0)
        axis.set_axisbelow(True)

    tick_positions = []
    tick_labels = []
    for index, chrom in enumerate(ordered_chroms):
        start = offsets[chrom]
        end = offsets[ordered_chroms[index + 1]] if index + 1 < len(ordered_chroms) else genome_length
        tick_positions.append((start + end) / 2.0)
        tick_labels.append(chrom.replace("chr", ""))
    baf_axis.set_xticks(tick_positions)
    baf_axis.set_xticklabels(tick_labels, fontsize=8)
    baf_axis.set_xlabel("Chromosome", fontsize=10)

    legend_entries = [
        Patch(facecolor=COLOUR_LOH, label="Region flagged LOH likely"),
        Patch(facecolor=COLOUR_EQUIVOCAL, label="Region equivocal"),
        Patch(facecolor=COLOUR_BALANCED, label="Region balanced"),
        Patch(facecolor=COLOUR_PANEL_SHADE, label="Panel window"),
        Patch(facecolor=COLOUR_PANEL_SHADE, hatch="///", edgecolor="#9AA5B1",
              label="Panel window, unassessable (<{0} het sites)".format(min_sites)),
    ]
    baf_axis.legend(handles=legend_entries, loc="lower left", fontsize=8,
                    ncol=5, frameon=True, framealpha=0.9)

    figure.text(
        0.5, 0.015,
        "BAF is measured only within panel windows. Unshaded genomic space carries no BAF "
        "information and must not be read as balanced heterozygosity.",
        ha="center", fontsize=8, color="#666666",
    )

    outpath = os.path.join(outdir, "{0}.genome_baf_cn.png".format(sample))
    figure.savefig(outpath, dpi=160, bbox_inches="tight")
    plt.close(figure)
    return outpath


# ---------------------------------------------------------------------------
# Figure 2: per-region detail
# ---------------------------------------------------------------------------

def plot_regions_detail(sample, screen_rows, regions, het_by_chrom,
                        cna_by_chrom, outdir, min_sites, max_panels=None):
    """
    Draw a grid of per-region BAF scatter plots with marginal histograms.

    Point-level display is what makes bimodality directly visible; the summary
    statistics in the screen table are a compression of exactly this picture,
    and reviewing them side by side is how a borderline call gets adjudicated.
    """
    selected = [region for region in regions if (sample, region["name"]) in screen_rows]
    if max_panels:
        # Prioritise flagged regions so a truncated figure still shows the signal.
        priority = {"LOH_LIKELY": 0, "EQUIVOCAL": 1, "NO_LOH": 2, "UNASSESSABLE": 3}
        selected.sort(key=lambda region: priority.get(
            screen_rows[(sample, region["name"])]["flag"], 4))
        selected = selected[:max_panels]

    if not selected:
        return None

    columns = 4
    rows_count = (len(selected) + columns - 1) // columns
    figure, axes = plt.subplots(
        rows_count, columns,
        figsize=(4.0 * columns, 2.7 * rows_count),
        squeeze=False,
    )

    for index, region in enumerate(selected):
        axis = axes[index // columns][index % columns]
        row = screen_rows[(sample, region["name"])]
        flag = row["flag"]
        colour = FLAG_COLOURS.get(flag, COLOUR_BALANCED)

        sites = [
            (position, minor)
            for position, _depth, minor in het_by_chrom.get(region["chrom"], [])
            if region["start"] <= position - 1 < region["end"]
        ]

        if flag == "UNASSESSABLE" or not sites:
            axis.set_facecolor("#F5F5F5")
            axis.text(
                0.5, 0.5,
                "Unassessable\n{0} het sites (minimum {1})".format(row["n_het"], min_sites),
                transform=axis.transAxes, ha="center", va="center",
                fontsize=9, color="#777777",
            )
            axis.set_xticks([])
            axis.set_yticks([])
        else:
            positions = [position / 1e6 for position, _minor in sites]
            values = [minor for _position, minor in sites]
            axis.scatter(positions, values, s=3.0, c=colour, alpha=0.5, linewidths=0)
            axis.axhline(0.5, color="#666666", linewidth=0.8)
            # Central band bounds mark where balanced heterozygous sites
            # are expected to sit; LOH shows as depletion between them.
            axis.axhspan(0.42, 0.58, color="#2E7D32", alpha=0.10, zorder=0)
            axis.set_ylim(0.0, 1.0)
            axis.tick_params(labelsize=7)
            axis.set_xlabel("Mb", fontsize=8)
            axis.set_ylabel("BAF", fontsize=8)

        subtitle_parts = ["n={0}".format(row["n_het"])]
        if row["median_dp"] is not None:
            subtitle_parts.append("dp={0:.0f}x".format(row["median_dp"]))
        if row["depletion_score"] is not None:
            subtitle_parts.append("CDR={0:.2f}".format(row["depletion_score"]))
        if row["bimodality"] is not None:
            subtitle_parts.append("bim={0:.1f}".format(row["bimodality"]))

        axis.set_title(
            "{0}\n{1}".format(region["name"], "  ".join(subtitle_parts)),
            fontsize=8.5, color=colour if flag != "NO_LOH" else "#333333",
        )

    # Blank any unused grid cells.
    for index in range(len(selected), rows_count * columns):
        axes[index // columns][index % columns].axis("off")

    figure.suptitle(
        "Sample {0} - per-region minor-allele BAF "
        "(green band: expected range for balanced heterozygosity)".format(sample),
        fontsize=12, y=0.997,
    )
    figure.tight_layout(rect=[0, 0, 1, 0.985])

    outpath = os.path.join(outdir, "{0}.region_baf.png".format(sample))
    figure.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(figure)
    return outpath


# ---------------------------------------------------------------------------
# Figure 3: cohort heatmap
# ---------------------------------------------------------------------------

def plot_cohort_heatmap(screen_rows, samples, regions, outdir, min_sites):
    """
    Draw a samples-by-regions heatmap of BAF deflection.

    Deflection is plotted rather than the categorical flag so that regions
    sitting just below the calling threshold remain visible; collapsing to the
    flag would hide exactly the borderline cases that most need review.
    Unassessable cells are hatched, never coloured, because any colour taken
    from the scale would be read as a measured value.
    """
    region_names = [region["name"] for region in regions
                    if any((sample, region["name"]) in screen_rows for sample in samples)]
    if not region_names or not samples:
        return None

    colormap = LinearSegmentedColormap.from_list(
        "deflection", ["#F7FBFF", "#FDD9A0", "#E8743B", "#8B1A0E"]
    )

    values = []
    for (_sample, _region), row in screen_rows.items():
        if row["flag"] != "UNASSESSABLE" and row["baf_deflection"] is not None:
            values.append(row["baf_deflection"])
    upper = max(0.10, statistics.median(values) + 4 * (statistics.pstdev(values) or 0.02)) if values else 0.30
    normaliser = Normalize(vmin=0.0, vmax=min(upper, 0.5))

    figure, axis = plt.subplots(
        figsize=(max(9.0, 0.30 * len(region_names)), max(3.0, 0.55 * len(samples) + 2.2))
    )

    for row_index, sample in enumerate(samples):
        for column_index, region_name in enumerate(region_names):
            row = screen_rows.get((sample, region_name))
            if row is None or row["flag"] == "UNASSESSABLE" or row["baf_deflection"] is None:
                axis.add_patch(plt.Rectangle(
                    (column_index, row_index), 1, 1,
                    facecolor="white", edgecolor="#BBBBBB",
                    hatch="///", linewidth=0.4,
                ))
                continue
            deflection = max(0.0, row["baf_deflection"])
            axis.add_patch(plt.Rectangle(
                (column_index, row_index), 1, 1,
                facecolor=colormap(normaliser(deflection)),
                edgecolor="#DDDDDD", linewidth=0.4,
            ))
            if row["flag"] == "LOH_LIKELY":
                # Mark called regions explicitly; colour alone should not carry
                # the categorical outcome.
                axis.plot(column_index + 0.5, row_index + 0.5, marker="o",
                          markersize=3.5, color="#111111")

    axis.set_xlim(0, len(region_names))
    axis.set_ylim(0, len(samples))
    axis.invert_yaxis()
    axis.set_xticks([index + 0.5 for index in range(len(region_names))])
    axis.set_xticklabels(region_names, rotation=90, fontsize=7)
    axis.set_yticks([index + 0.5 for index in range(len(samples))])
    axis.set_yticklabels(samples, fontsize=9)
    axis.set_title(
        "Cohort BAF deflection by panel region "
        "(dot: flagged LOH likely; hatched: unassessable, <{0} het sites)".format(min_sites),
        fontsize=11, pad=12,
    )

    mappable = plt.cm.ScalarMappable(cmap=colormap, norm=normaliser)
    mappable.set_array([])
    colourbar = figure.colorbar(mappable, ax=axis, fraction=0.025, pad=0.015)
    colourbar.set_label("Mean |BAF - 0.5| (allelic imbalance)", fontsize=9)
    colourbar.ax.tick_params(labelsize=8)

    figure.tight_layout()
    outpath = os.path.join(outdir, "cohort_baf_deflection_heatmap.png")
    figure.savefig(outpath, dpi=160, bbox_inches="tight")
    plt.close(figure)
    return outpath


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def read_plot_inputs(path):
    """Read the plotting input table: sample, vcf_path, cna_seg, params_txt."""
    entries = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        index = {name.strip(): position for position, name in enumerate(header)}
        required = ["sample", "vcf_path"]
        for column in required:
            if column not in index:
                raise ValueError("Plot input table missing required column: {0}".format(column))
        for line in handle:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split("\t")

            def value(name):
                if name not in index or index[name] >= len(fields):
                    return None
                raw = fields[index[name]].strip()
                return None if raw in ("", "NA") else raw

            entries.append({
                "sample": value("sample"),
                "vcf_path": value("vcf_path"),
                "cna_seg": value("cna_seg"),
                "params_txt": value("params_txt"),
            })
    return entries


def build_argument_parser():
    parser = argparse.ArgumentParser(
        description="Combined BAF and copy number visualisation for panel adaptive sampling data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--screen", required=True,
                        help="Output TSV from baf_loh_screen.py.")
    parser.add_argument("--bed", required=True, help="Panel BED file.")
    parser.add_argument("--sample-map", dest="sample_map", required=True,
                        help="Table of sample, vcf_path, cna_seg, params_txt.")
    parser.add_argument("--outdir", required=True, help="Output directory for figures.")
    parser.add_argument("--min-site-depth", dest="min_site_depth", type=int, default=8,
                        help="Minimum het site depth, matching the screen settings.")
    parser.add_argument("--min-sites", dest="min_sites", type=int, default=30,
                        help="Assessability threshold, matching the screen settings.")
    parser.add_argument("--max-region-panels", dest="max_region_panels", type=int, default=None,
                        help="Cap on per-region detail panels; flagged regions are kept first.")
    parser.add_argument("--skip-genome", action="store_true",
                        help="Skip the genome-wide figure.")
    parser.add_argument("--skip-regions", action="store_true",
                        help="Skip the per-region detail figure.")
    return parser


def main(argv=None):
    parser = build_argument_parser()
    args = parser.parse_args(argv)

    os.makedirs(args.outdir, exist_ok=True)

    regions = read_bed(args.bed)
    screen_rows, screen_samples = read_screen_table(args.screen)
    inputs = read_plot_inputs(args.sample_map)

    sys.stderr.write("Loaded {0} regions, {1} screened samples, {2} plot inputs\n".format(
        len(regions), len(screen_samples), len(inputs)))

    written = []
    plotted_samples = []

    for entry in inputs:
        sample = entry["sample"]
        if sample is None:
            continue
        plotted_samples.append(sample)

        vcf_paths = collect_vcf_paths(entry["vcf_path"]) if entry["vcf_path"] else []
        het_by_chrom = load_het_sites(vcf_paths, args.min_site_depth) if vcf_paths else {}
        cna_by_chrom = read_cna_seg(entry["cna_seg"])
        params = read_ichor_params(entry["params_txt"])

        if not cna_by_chrom:
            sys.stderr.write(
                "{0}: no copy number bins loaded; BAF panel will be drawn without a CN track. "
                "Copy-neutral LOH cannot be distinguished from hemizygous loss without it.\n".format(sample)
            )

        if not args.skip_genome:
            path = plot_genome_wide(sample, screen_rows, regions, het_by_chrom,
                                    cna_by_chrom, params, args.outdir, args.min_sites)
            if path:
                written.append(path)

        if not args.skip_regions:
            path = plot_regions_detail(sample, screen_rows, regions, het_by_chrom,
                                       cna_by_chrom, args.outdir, args.min_sites,
                                       args.max_region_panels)
            if path:
                written.append(path)

    heatmap_samples = plotted_samples or screen_samples
    path = plot_cohort_heatmap(screen_rows, heatmap_samples, regions,
                               args.outdir, args.min_sites)
    if path:
        written.append(path)

    for path in written:
        sys.stderr.write("Wrote {0}\n".format(path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
