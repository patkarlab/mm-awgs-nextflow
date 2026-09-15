"""
parsers/copy_number.py

Depth-based copy number from ichorCNA, by way of ichorkaryo.

Distinct from parsers/ichor.py, which reads ichorCNA's own output directly for
the tumour-fraction panel. This reads the karyotype JSON ichorkaryo derives
from the same segments: arm calls, an ISCN string, a cytoband table, per-panel
gene rows, and the QC that says whether a negative means anything.

What is and is not here
-----------------------
Allele-specific copy number is not. ichorCNA is depth-only and cannot separate
2:0 from 1:1, and allelic state is now measured at the eight wide panel windows
by genebaf rather than genome-wide. Minor allele copy number and copy-neutral
LOH therefore do not appear in this tab.

Purity comes from flow cytometry or not at all. ichorCNA's own tumour fraction
is the magnitude of copy-number deviation, not the proportion of tumour present
in sorted plasma cells, and it is surfaced under its own name as CNA burden so
it cannot be read as purity.

Reference frames
----------------
Two are in play and they are not the same. Gains and losses are called against
the genome's modal copy number, so on a near-triploid sample three copies is
neutral. Trisomies and the ploidy class are referenced to a constitutional two,
because hyperdiploidy is defined against a normal karyotype. ichorkaryo emits a
warning whenever the two differ and it is passed through.

Expected layout, as assembled by build_report_bundle.sh:

    <bundle>/<sample>/copy_number/<sample>.ichorkaryo.json
                                  <sample>.cn_genome.png
                                  <sample>.cn_grid.png
                                  <sample>.cn_chr<N>.png
"""

from __future__ import annotations

import json
import os

SUBDIR = "copy_number"

# The published name is <sample>.ichorkaryo.json. The others are accepted so a
# bundle assembled before the rename, or by hand, still renders.
KARYOTYPE_NAMES = ("%s.ichorkaryo.json", "%s.karyotype.json",
                   "ichorkaryo.json", "karyotype.json")

CHROM_ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]


def _first(directory, *names):
    for name in names:
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            return candidate
    return None


def _arm_sort_key(label):
    """Karyotype order for arm labels, so 10p does not sort before 2p."""
    body = label.rstrip("pq")
    suffix = label[len(body):]
    index = {"X": 23, "Y": 24}.get(body)
    if index is None:
        try:
            index = int(body)
        except ValueError:
            index = 99
    return (index, {"p": 0, "q": 1}.get(suffix, 2))


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse(effective_dir, sample):
    directory = os.path.join(effective_dir, SUBDIR)
    result = {
        "found": False,
        "searched": directory,
        "headline": {},
        "figures": {},
        "qc": {},
        "arms": [],
        "cytobands": [],
        "altered_segments": [],
        "genes": [],
        "chromosomes": [],
        "warnings": [],
    }

    candidates = [n % sample if "%s" in n else n for n in KARYOTYPE_NAMES]
    karyotype_path = _first(directory, *candidates)
    if not karyotype_path:
        return result

    try:
        with open(karyotype_path) as handle:
            data = json.load(handle)
    except (ValueError, OSError):
        result["warnings"].append("could not parse %s" % karyotype_path)
        return result

    result["found"] = True
    fit = data.get("fit", {}) or {}
    qc = data.get("qc", {}) or {}

    purity = _as_float(fit.get("purity"))
    purity_supplied = purity is not None and fit.get("purity_source") == "flow"

    result["headline"] = {
        "iscn": data.get("ISCN_karyotype", ""),
        "sex": data.get("sex", ""),
        "purity": purity,
        "purity_source": ("flow cytometry" if purity_supplied
                          else (fit.get("purity_source") or "not supplied")),
        "purity_supplied": purity_supplied,
        "ploidy": _as_float(fit.get("ploidy")),
        "ploidy_class": data.get("ploidy_class"),
        "hyperdiploid": bool(data.get("hyperdiploid")),
        "trisomies": data.get("trisomies") or [],
        "canonical_trisomies": data.get("canonical_trisomies") or [],
        "chromosome_count": data.get("estimated_chromosome_count"),
        "baseline_cn": fit.get("baseline_copy_number"),
        # ichorCNA's fitted fraction, named for what it measures on sorted
        # plasma cells. Never a purity.
        "cna_burden": _as_float(fit.get("cna_burden")),
        "n_segments": data.get("n_segments"),
    }

    # A detection limit above the FISH percentages being compared against means
    # a negative call is not yet meaningful, so it belongs in the headline QC
    # rather than buried.
    bins_total = qc.get("n_bins")
    bins_usable = qc.get("n_bins_usable")
    excluded = None
    if bins_total:
        excluded = 1.0 - (bins_usable or 0) / float(bins_total)

    result["qc"] = {
        "detection_limit_5mb": _as_float(qc.get("detection_limit_5mb")),
        "median_reads_per_bin": qc.get("median_reads_per_bin"),
        # Corrected is the one near 1 when counting-limited; raw is before
        # ichorCNA's GC and mappability correction and is not comparable to
        # the figure mmkaryo used to report.
        "overdispersion": _as_float(qc.get("overdispersion")),
        "overdispersion_raw": _as_float(qc.get("overdispersion_raw_counts")),
        "sigma_log2_per_bin": _as_float(qc.get("sigma_log2_per_bin")),
        "bin_size": qc.get("bin_size"),
        "bins_total": bins_total,
        "bins_usable": bins_usable,
        "genome_excluded_fraction": excluded,
        "baseline_method": fit.get("baseline_source"),
    }

    # ichorkaryo already decides what is worth warning about, including the
    # mixed reference frames and a baseline covering too little of the genome.
    result["warnings"] = list(data.get("warnings") or [])
    if not purity_supplied:
        result["warnings"].append(
            "No flow cytometry purity was supplied, so cell fractions are not "
            "reported. The CNA burden shown is ichorCNA's fitted deviation, "
            "not the proportion of tumour present.")
    if result["qc"]["detection_limit_5mb"] is None:
        result["warnings"].append(
            "No detection limit could be computed, so a negative call here "
            "cannot be interpreted.")

    # Arms arrive keyed by label. The template wants rows, in karyotype order.
    arms = data.get("arms") or {}
    result["arms"] = [{
        "arm": label,
        "chromosome": (arms[label] or {}).get("chrom", ""),
        "log2_ratio": (arms[label] or {}).get("log2"),
        "copy_number": (arms[label] or {}).get("cn"),
        "event": (arms[label] or {}).get("event", ""),
        # Fraction of the arm the majority call covers. A partial event reads
        # as partial rather than being rounded into a whole-arm statement.
        "arm_fraction": (arms[label] or {}).get("arm_fraction"),
        "cell_fraction": (arms[label] or {}).get("cell_fraction"),
    } for label in sorted(arms, key=_arm_sort_key)]

    # Already filtered to altered bands by ichorkaryo, and already carrying the
    # fields the table renders.
    result["cytobands"] = data.get("cytobands") or []
    result["altered_segments"] = data.get("altered_segments") or []
    result["genes"] = data.get("genes") or []

    for name, key in (("cn_genome.png", "genome"), ("cn_grid.png", "grid")):
        path = _first(directory, "%s.%s" % (sample, name), name)
        if path:
            result["figures"][key] = os.path.join(SUBDIR,
                                                  os.path.basename(path))

    # Per-chromosome pages, in karyotype order rather than the lexical order a
    # directory listing gives, which would put chr10 before chr2.
    for chrom in CHROM_ORDER:
        path = _first(directory, "%s.cn_%s.png" % (sample, chrom),
                      "cn_%s.png" % chrom)
        if path:
            result["chromosomes"].append({
                "chrom": chrom,
                "label": chrom[3:],
                "src": os.path.join(SUBDIR, os.path.basename(path)),
            })

    return result
