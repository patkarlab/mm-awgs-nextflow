#!/usr/bin/env bash
#
# build_report_bundle.sh
# ----------------------
# Build a per-sample, downloadable reporting bundle from a pipeline results
# directory. Copies REAL files (not symlinks) so the resulting tarball is
# self-contained and portable off gandalf.
#
# Layout produced:
#   <bundle>/
#     <SAMPLE>/
#       snv/            <SAMPLE>.clinical.tsv, <SAMPLE>.filtered.tsv
#       translocations/ <SAMPLE>.mm_annotated.tsv, <SAMPLE>.translocations.tsv
#       cnv/            <SAMPLE>.ichor_all_sols.pdf, <SAMPLE>.ichor_params.txt
#     v6_filter_summary.tsv        (cohort-level, if present)
#   <bundle>.zip
#
# The ".v6_" label from the pipeline is stripped in the bundle copies only;
# source files are untouched. (The panel is v7; the v6 filename is legacy.)
#
# Usage:
#   ./build_report_bundle.sh <results_dir> [bundle_name]
# Example:
#   ./build_report_bundle.sh results_v7_20260713_24h report_v7_20260713_24h
#
set -euo pipefail

RESULTS="${1:?Usage: build_report_bundle.sh <results_dir> [bundle_name]}"
BUNDLE="${2:-report_$(basename "$RESULTS")}"

if [[ ! -d "$RESULTS" ]]; then
  echo "ERROR: results dir not found: $RESULTS" >&2
  exit 1
fi

# Derive the sample list from the annotated TSVs (most reliable anchor).
mapfile -t SAMPLES < <(
  find "$RESULTS" -name '*.mm_annotated.tsv' -printf '%f\n' 2>/dev/null \
    | sed 's/\.mm_annotated\.tsv$//' | sort -u
)
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  # fallback: derive from clinical TSVs
  mapfile -t SAMPLES < <(
    find "$RESULTS" -name '*clinical.tsv' -printf '%f\n' 2>/dev/null \
      | sed -E 's/(\.somatic_candidates)?\.v6_clinical\.tsv$//' | sort -u
  )
fi
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no samples found (no *.mm_annotated.tsv or *clinical.tsv under $RESULTS)" >&2
  exit 1
fi

echo "Samples detected: ${SAMPLES[*]}"
rm -rf "$BUNDLE"
mkdir -p "$BUNDLE"

# copy_first <destdir> <destname> <find-pattern...>
# Finds the first matching file under RESULTS filtered to the sample, copies it.
copy_first() {
  local destdir="$1"; shift
  local destname="$1"; shift
  local sample="$1"; shift
  local src
  src=$(find "$RESULTS" -type f "$@" 2>/dev/null | grep -F "$sample" | head -1 || true)
  if [[ -n "$src" && -f "$src" ]]; then
    mkdir -p "$destdir"
    cp -L "$src" "$destdir/$destname"
    echo "  + $destname"
  else
    echo "  - (missing) $destname" >&2
  fi
}

for s in "${SAMPLES[@]}"; do
  echo "== $s =="
  d="$BUNDLE/$s"

  # SNV (strip the v6 label in the bundle)
  copy_first "$d/snv" "${s}.clinical.tsv" "$s" -name '*clinical.tsv'
  copy_first "$d/snv" "${s}.filtered.tsv" "$s" -name '*filtered.tsv'

  # Translocations
  copy_first "$d/translocations" "${s}.mm_annotated.tsv"   "$s" -name '*.mm_annotated.tsv'
  copy_first "$d/translocations" "${s}.translocations.tsv" "$s" -name '*.translocations.tsv'

  # CNV: ONLY the all_sols PDF, plus the params (tumor fraction / ploidy).
  copy_first "$d/cnv" "${s}.ichor_all_sols.pdf" "$s" -path '*ichor*' -name '*all_sols*.pdf'
  copy_first "$d/cnv" "${s}.ichor_params.txt"   "$s" -path '*ichor*' -name '*params.txt'

  # QC: on-target panel-region coverage (table + chart) and read-length/qscore summary.
  copy_first "$d/qc" "${s}.region_coverage.tsv" "$s" -name '*region_coverage.tsv'
  copy_first "$d/qc" "${s}.region_coverage.png" "$s" -name '*region_coverage.png'
  copy_first "$d/qc" "${s}.readlen_qscore.tsv"  "$s" -name '*readlen_qscore.tsv'
  copy_first "$d/qc" "${s}.readlen_hist.png"    "$s" -name '*readlen_hist.png'
  copy_first "$d/qc" "${s}.qscore_hist.png"     "$s" -name '*qscore_hist.png'
done

# Cohort-level summary (single file, not per-sample).
summary=$(find "$RESULTS" -name 'v6_filter_summary.tsv' 2>/dev/null | head -1 || true)
if [[ -n "$summary" ]]; then
  cp -L "$summary" "$BUNDLE/filter_summary.tsv"
  echo "+ cohort filter_summary.tsv"
fi

# Tar it up.
if command -v zip >/dev/null 2>&1; then
  zip -rq "${BUNDLE}.zip" "$BUNDLE"
else
  echo "WARNING: 'zip' not found; falling back to tar.gz" >&2
  tar czf "${BUNDLE}.tar.gz" "$BUNDLE"
fi
echo ""
echo "Bundle tree: $BUNDLE/"
if [ -f "${BUNDLE}.zip" ]; then echo "Archive:     ${BUNDLE}.zip"; else echo "Archive:     ${BUNDLE}.tar.gz"; fi
du -sh "${BUNDLE}.zip" 2>/dev/null || du -sh "${BUNDLE}.tar.gz"
echo ""
echo "Contents:"
find "$BUNDLE" -type f | sed "s|^$BUNDLE/||" | sort
