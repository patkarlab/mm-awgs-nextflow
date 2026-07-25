#!/usr/bin/env bash
#
# make_igv_snapshots.sh
#
# Generate IGV snapshots for a results directory that has already been
# produced by the pipeline. This exists so IGV can be added to runs that
# completed before the IGV stage was wired into the workflow, without
# re-running any calling.
#
# Two evidence classes are rendered per sample:
#   somatic         one page covering the clinical SNV table, against hg38
#   translocations  two standalone pages per event (both breakpoints),
#                   against T2T, plus a manifest the dashboard reads
#
# Usage:
#   tools/make_igv_snapshots.sh <results_dir> [flanking]
#
# Example:
#   tools/make_igv_snapshots.sh results_v7_20260713_24h
#   tools/make_igv_snapshots.sh results_v7_20260713_24h 8000
#
# Overrides, if reference paths differ from the defaults below:
#   T2T_REF=/path/chm13v2.0.ucsc.fa HG38_REF=/path/hg38.fasta \
#     tools/make_igv_snapshots.sh results_v7_20260713_24h
#
# Set DRY_RUN=1 to print the resolved inputs per sample and exit without
# rendering anything. Run that first on a new results layout.

set -euo pipefail

# nohup and Nextflow do not inherit an interactive shell environment, so the
# environment is activated explicitly.
source /home/hemat/anaconda3/etc/profile.d/conda.sh
conda activate awgs_sv

RESULTS="${1:?Usage: make_igv_snapshots.sh <results_dir> [flanking]}"
FLANKING="${2:-5000}"
DRY_RUN="${DRY_RUN:-0}"

T2T_REF="${T2T_REF:-/goast/nikhil_awgs_testing/t2t/refs/chm13v2.0.ucsc.fa}"
HG38_REF="${HG38_REF:-/goast/hemat_data/references/hg38_broad/Homo_sapiens_assembly38.fasta}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IGV_SNAPSHOTS="${IGV_SNAPSHOTS:-${SCRIPT_DIR}/../bin/igv_snapshots.py}"

if [[ ! -d "$RESULTS" ]]; then
  echo "ERROR: results directory not found: $RESULTS" >&2
  exit 1
fi
if [[ ! -f "$IGV_SNAPSHOTS" ]]; then
  echo "ERROR: igv_snapshots.py not found at $IGV_SNAPSHOTS" >&2
  exit 1
fi
for ref in "$T2T_REF" "$HG38_REF"; do
  if [[ ! -f "$ref" ]]; then
    echo "ERROR: reference not found: $ref" >&2
    echo "       override with T2T_REF= / HG38_REF=" >&2
    exit 1
  fi
done

# Samples are derived from the annotated translocation tables, which are the
# most reliable per-sample anchor in this tree.
mapfile -t SAMPLES < <(
  find "$RESULTS" -name '*.mm_annotated.tsv' -printf '%f\n' 2>/dev/null \
    | sed 's/\.mm_annotated\.tsv$//' | sort -u
)
if [[ ${#SAMPLES[@]} -eq 0 ]]; then
  echo "ERROR: no *.mm_annotated.tsv found under $RESULTS" >&2
  exit 1
fi

echo "Samples:  ${SAMPLES[*]}"
echo "Flanking: ${FLANKING} bp"
echo "T2T ref:  ${T2T_REF}"
echo "hg38 ref: ${HG38_REF}"
echo ""

# find_one <description> <find-args...>
# Prints the first match, or nothing. Never fails the script; a missing input
# is reported per sample and that evidence class is skipped.
find_one() {
  find "$RESULTS" -type f "$@" 2>/dev/null | head -1 || true
}

rc=0
for sample in "${SAMPLES[@]}"; do
  echo "== ${sample} =="
  outdir="${RESULTS}/igv/${sample}"

  mm_tsv=$(find_one -name "${sample}.mm_annotated.tsv")
  clin_tsv=$(find_one -name "${sample}*clinical.tsv" -size +1c)

  # BAMs: prefer the reference-specific subtree, fall back to any BAM whose
  # name carries the sample identifier.
  t2t_bam=$(find "$RESULTS" -type f -ipath '*t2t*' -name "${sample}*.bam" 2>/dev/null | head -1 || true)
  [[ -z "$t2t_bam" ]] && t2t_bam=$(find_one -name "${sample}*t2t*.bam")
  hg38_bam=$(find "$RESULTS" -type f -ipath '*hg38*' -name "${sample}*.bam" 2>/dev/null | head -1 || true)
  [[ -z "$hg38_bam" ]] && hg38_bam=$(find_one -name "${sample}*hg38*.bam")

  echo "  translocations tsv : ${mm_tsv:-MISSING}"
  echo "  clinical tsv       : ${clin_tsv:-MISSING or empty}"
  echo "  T2T bam            : ${t2t_bam:-MISSING}"
  echo "  hg38 bam           : ${hg38_bam:-MISSING}"

  if [[ "$DRY_RUN" == "1" ]]; then
    echo ""
    continue
  fi

  mkdir -p "${outdir}/translocations" "${outdir}/somatic"

  if [[ -n "$mm_tsv" && -n "$t2t_bam" ]]; then
    echo "  rendering translocation breakpoints..."
    if ! python3 "$IGV_SNAPSHOTS" \
        --mode translocations \
        --sample "$sample" \
        --sites-tsv "$mm_tsv" \
        --bam "$t2t_bam" \
        --fasta "$T2T_REF" \
        --out-html "${outdir}/translocations/${sample}.translocations.html" \
        --out-dir "${outdir}/translocations" \
        --flanking "$FLANKING"; then
      echo "  WARNING: translocation snapshots failed for ${sample}" >&2
      rc=1
    fi
  else
    echo "  skipping translocations (missing table or T2T BAM)"
  fi

  if [[ -n "$clin_tsv" && -n "$hg38_bam" ]]; then
    echo "  rendering clinical SNVs..."
    if ! python3 "$IGV_SNAPSHOTS" \
        --mode somatic \
        --sample "$sample" \
        --sites-tsv "$clin_tsv" \
        --bam "$hg38_bam" \
        --fasta "$HG38_REF" \
        --out-html "${outdir}/somatic/${sample}.somatic.html" \
        --flanking "$FLANKING"; then
      echo "  WARNING: somatic snapshots failed for ${sample}" >&2
      rc=1
    fi
  else
    echo "  skipping somatic (missing clinical table or hg38 BAM)"
  fi

  echo ""
done

echo "Output tree: ${RESULTS}/igv/"
find "${RESULTS}/igv" -name '*.html' 2>/dev/null | wc -l | xargs echo "HTML pages:"
du -sh "${RESULTS}/igv" 2>/dev/null || true

exit "$rc"
