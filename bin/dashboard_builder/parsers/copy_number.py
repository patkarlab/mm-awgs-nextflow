"""
parsers/copy_number.py

Allele-specific copy number from the T2T track (mmkaryo / mmbaf / mmplot).

Distinct from parsers/ichor.py. ichorCNA is a tumour-fraction estimator that
emits copy number from off-target reads against hg38; this reads segmented,
allele-specific copy number against T2T, with a clonal cell fraction per
segment and an explicit detection limit.

Two fields are deliberately NOT surfaced, because both are documented as
unreliable until BAF gets its own segmentation and would read as authoritative
in a report:

  - the baseline disomic/trisomic verdict from mmbaf, which is confounded by
    theta inheriting depth segmentation
  - any purity fitted with purity left free, which is not identifiable from
    read depth alone

Where purity came from IS reported, so a reader can see whether the copy
numbers rest on a flow measurement or on a depth fit.

Expected layout, as assembled by build_report_bundle.sh:

    <bundle>/<sample>/copy_number/<sample>.cn_genome.png
                                  <sample>.cn_grid.png
                                  <sample>.karyotype.json
                                  <sample>.cytobands.tsv
                                  <sample>.arms.tsv
                                  <sample>.segments.baf.tsv
"""

from __future__ import annotations

import json
import os

SUBDIR = "copy_number"


def _read_tsv(path, limit=None):
    if not path or not os.path.exists(path):
        return [], []
    rows = []
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < len(header):
                continue
            rows.append(dict(zip(header, fields)))
            if limit and len(rows) >= limit:
                break
    return header, rows


def _first(directory, *names):
    for name in names:
        candidate = os.path.join(directory, name)
        if os.path.exists(candidate):
            return candidate
    return None


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
        "chromosomes": [],
        "warnings": [],
    }

    karyotype_path = _first(directory, "%s.karyotype.json" % sample,
                            "karyotype.json")
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
    hyperdiploidy = data.get("hyperdiploidy", {}) or {}

    purity_supplied = bool(fit.get("purity_supplied"))
    result["headline"] = {
        "iscn": data.get("ISCN_karyotype", ""),
        "sex": data.get("inferred_sex", ""),
        "purity": fit.get("purity"),
        "purity_source": "flow cytometry" if purity_supplied
                         else "fitted from read depth",
        "purity_supplied": purity_supplied,
        "ploidy": fit.get("ploidy"),
        "hyperdiploid": hyperdiploidy.get("called"),
        "trisomies": hyperdiploidy.get("trisomic_chromosomes") or [],
        "n_segments": data.get("n_segments"),
    }

    # A detection limit above the FISH percentages being compared against means
    # a negative call is not yet meaningful, so it belongs in the headline QC
    # rather than buried.
    result["qc"] = {
        "detection_limit_5mb": qc.get("detection_limit_5mb"),
        "median_reads_per_bin": qc.get("median_reads_per_bin"),
        "overdispersion": qc.get("overdispersion"),
        "bins_retained": qc.get("bins_retained"),
        "genome_excluded_fraction": qc.get("genome_excluded_fraction"),
        "baseline_method": qc.get("baseline_method"),
    }

    if not purity_supplied:
        result["warnings"].append(
            "Purity was fitted from read depth, which cannot separate purity "
            "from ploidy. Supply the post-sort flow value in the sample sheet "
            "to pin it.")
    if fit.get("poor_fit"):
        result["warnings"].append(
            "Segments are not resolving to integer copy states; treat the "
            "integer calls as provisional and read the cell fractions.")
    unexplained = _as_float(fit.get("unexplained_fraction"))
    if unexplained is not None and unexplained > 0.10:
        result["warnings"].append(
            "%.0f%% of the genome is not explained by a single clone. Per-"
            "segment cell fractions are the meaningful readout there, not the "
            "integer copy numbers." % (100 * unexplained))

    result["altered_segments"] = data.get("altered_segments", []) or []

    for name, key in (("cn_genome.png", "genome"), ("cn_grid.png", "grid")):
        path = _first(directory, "%s.%s" % (sample, name), name)
        if path:
            result["figures"][key] = os.path.join(SUBDIR,
                                                  os.path.basename(path))

    # Per-chromosome pages, in karyotype order rather than the lexical order
    # a directory listing gives (which would put chr10 before chr2).
    order = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]
    for chrom in order:
        path = _first(directory, "%s.cn_%s.png" % (sample, chrom),
                      "cn_%s.png" % chrom)
        if path:
            result["chromosomes"].append({
                "chrom": chrom,
                "label": chrom[3:],
                "src": os.path.join(SUBDIR, os.path.basename(path)),
            })

    _header, arms = _read_tsv(_first(directory, "%s.arms.tsv" % sample))
    result["arms"] = arms

    _header, bands = _read_tsv(_first(directory, "%s.cytobands.tsv" % sample))
    # Only altered bands are worth tabulating; a full band list is thousands of
    # rows of copy number 2.
    result["cytobands"] = [b for b in bands
                           if b.get("copy_number") not in (None, "", "2")]

    return result
