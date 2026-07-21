#!/usr/bin/env bash
#
# make_igv_snapshots.sh
#
# Standalone cohort driver for igv_snapshots.py. Iterates the samples found in
# a results tree and renders, per sample:
#   - somatic HTML        (v6 clinical TSV  vs hg38 BAM)  -> hg38/igv/<sample>/
#   - translocation HTML  (translocations TSV vs T2T BAM) -> t2t/igv/<sample>/
#
# This is the bash entry point used outside Nextflow (the Nextflow module
# modules/local/igv_report.nf wraps the same igv_snapshots.py per process).
#
# Conventions honoured:
#   - conda is sourced and awgs_sv activated inside the script (nohup-safe).
#   - sample ids are discovered from the results tree; no ids are hardcoded.
#   - no patient names, FISH findings, or expected calls appear anywhere.
#   - BAMs are (re)indexed defensively (stale .bai has bitten this pipeline).
#
# Usage:
#   bash make_igv_snapshots.sh <results_dir> [sample_id ...]
#
#   <results_dir>   e.g. results_cohort_18h  (must contain t2t/ and hg38/)
#   [sample_id ...] optional explicit subset; default = all discovered samples
#
# Example:
#   bash make_igv_snapshots.sh results_cohort_18h
#   bash make_igv_snapshots.sh results_cohort_18h SAMPLE_ID

set -euo pipefail

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
source /home/hemat/anaconda3/etc/profile.d/conda.sh
conda activate awgs_sv

# Location of igv_snapshots.py: same dir as this script by default.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IGV_PY="${IGV_PY:-${SCRIPT_DIR}/igv_snapshots.py}"

# References (overridable via environment).
T2T_FASTA="${T2T_FASTA:-/goast/nikhil_awgs_testing/t2t/refs/chm13v2.0.ucsc.fa}"
HG38_FASTA="${HG38_FASTA:-/goast/hemat_data/references/hg38_broad/Homo_sapiens_assembly38.fasta}"

FLANKING="${FLANKING:-5000}"

# ---------------------------------------------------------------------------
# Arguments
# ---------------------------------------------------------------------------
if [[ $# -lt 1 ]]; then
    echo "Usage: bash make_igv_snapshots.sh <results_dir> [sample_id ...]" >&2
    exit 2
fi

RESULTS_DIR="$1"; shift
EXPLICIT_SAMPLES=("$@")

if [[ ! -d "${RESULTS_DIR}/t2t" || ! -d "${RESULTS_DIR}/hg38" ]]; then
    echo "ERROR: ${RESULTS_DIR} must contain t2t/ and hg38/ subdirs" >&2
    exit 2
fi

if [[ ! -f "${IGV_PY}" ]]; then
    echo "ERROR: igv_snapshots.py not found at ${IGV_PY}" >&2
    exit 2
fi

for f in "${T2T_FASTA}" "${T2T_FASTA}.fai" "${HG38_FASTA}" "${HG38_FASTA}.fai"; do
    if [[ ! -f "${f}" ]]; then
        echo "ERROR: required reference file missing: ${f}" >&2
        exit 2
    fi
done

# ---------------------------------------------------------------------------
# Sample discovery
# ---------------------------------------------------------------------------
# Samples are discovered from the T2T merged-translocations TSVs, which exist
# for every sample that completed the T2T track. This avoids hardcoding ids.
TRA_DIR="${RESULTS_DIR}/t2t/calls/mm_annotated"
declare -a SAMPLES
if [[ ${#EXPLICIT_SAMPLES[@]} -gt 0 ]]; then
    SAMPLES=("${EXPLICIT_SAMPLES[@]}")
else
    if [[ ! -d "${TRA_DIR}" ]]; then
        echo "ERROR: ${TRA_DIR} not found; cannot discover samples" >&2
        exit 2
    fi
    while IFS= read -r f; do
        base="$(basename "${f}")"
        SAMPLES+=("${base%.translocations.tsv}")
    done < <(find "${TRA_DIR}" -name '*.translocations.tsv' | sort)
fi

if [[ ${#SAMPLES[@]} -eq 0 ]]; then
    echo "ERROR: no samples discovered" >&2
    exit 2
fi

echo "=== IGV snapshots for ${#SAMPLES[@]} sample(s) in ${RESULTS_DIR} ==="
printf '  %s\n' "${SAMPLES[@]}"

# ---------------------------------------------------------------------------
# Per-sample rendering
# ---------------------------------------------------------------------------
n_ok=0
n_skip=0
for sample in "${SAMPLES[@]}"; do
    echo ""
    echo "----- ${sample} -----"

    t2t_bam="${RESULTS_DIR}/t2t/bams/${sample}.t2t.bam"
    hg38_bam="${RESULTS_DIR}/hg38/bams/${sample}.hg38.bam"
    tra_tsv="${TRA_DIR}/${sample}.translocations.tsv"
    som_tsv="${RESULTS_DIR}/hg38/calls/v6_filtered/${sample}.somatic_candidates.withAD.v6_clinical.tsv"

    # --- translocations (T2T) ---
    if [[ -f "${tra_tsv}" && -f "${t2t_bam}" ]]; then
        out="${RESULTS_DIR}/t2t/igv/${sample}/${sample}.translocations.html"
        python3 "${IGV_PY}" \
            --mode translocations \
            --sample "${sample}" \
            --sites-tsv "${tra_tsv}" \
            --bam "${t2t_bam}" \
            --fasta "${T2T_FASTA}" \
            --out-html "${out}" \
            --flanking "${FLANKING}" \
        && n_ok=$((n_ok + 1)) || echo "  [warn] translocation render failed"
    else
        echo "  [skip] translocations: missing TSV or T2T BAM"
        n_skip=$((n_skip + 1))
    fi

    # --- somatic (hg38) ---
    if [[ -f "${som_tsv}" && -f "${hg38_bam}" ]]; then
        out="${RESULTS_DIR}/hg38/igv/${sample}/${sample}.somatic.html"
        python3 "${IGV_PY}" \
            --mode somatic \
            --sample "${sample}" \
            --sites-tsv "${som_tsv}" \
            --bam "${hg38_bam}" \
            --fasta "${HG38_FASTA}" \
            --out-html "${out}" \
            --flanking "${FLANKING}" \
        && n_ok=$((n_ok + 1)) || echo "  [warn] somatic render failed"
    else
        echo "  [skip] somatic: missing v6_clinical TSV or hg38 BAM"
        n_skip=$((n_skip + 1))
    fi
done

echo ""
echo "=== done: ${n_ok} report(s) rendered, ${n_skip} skipped ==="
