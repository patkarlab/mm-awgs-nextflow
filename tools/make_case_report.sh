#!/usr/bin/env bash
# make_case_report.sh
#
# Assemble a self-contained per-case report folder from an mm-awgs-nextflow
# results tree. For each sample it gathers:
#   - ichorCNA genome-wide PDFs (incl. *_genomeWide_all_sols.pdf) + .params.txt
#   - the v6 somatic report we generated (.v6_clinical.tsv + .v6_filtered.tsv)
#   - the T2T SV summary (.mm_annotated.tsv)
#
# Files are COPIED (not symlinked) so the report folder is portable.
#
# Usage:
#   bash make_case_report.sh <results_dir> [report_dir]
#
#   <results_dir>  pipeline results dir (the one containing hg38/ and t2t/)
#   [report_dir]   output; defaults to <results_dir>/report
#
# Layout produced:
#   report/<sample>/ichorcna/   *_genomeWide*.pdf, *.params.txt
#   report/<sample>/somatic/    <sample>.v6_clinical.tsv, <sample>.v6_filtered.tsv
#   report/<sample>/sv/         <sample>.mm_annotated.tsv
#   report/MANIFEST.tsv         per-sample presence check

set -euo pipefail

RESULTS_DIR="${1:-}"
if [[ -z "${RESULTS_DIR}" ]]; then
    echo "Usage: bash make_case_report.sh <results_dir> [report_dir]" >&2
    exit 1
fi
REPORT_DIR="${2:-${RESULTS_DIR}/report}"

ICHOR_BASE="${RESULTS_DIR}/hg38/calls/ichorcna"
V6_BASE="${RESULTS_DIR}/hg38/calls/v6_filtered"
MM_BASE="${RESULTS_DIR}/t2t/calls/mm_annotated"

# Discover samples from the v6 clinical reports (our anchor artifact).
shopt -s nullglob
CLIN=( "${V6_BASE}"/*.v6_clinical.tsv )
shopt -u nullglob
if [[ ${#CLIN[@]} -eq 0 ]]; then
    echo "ERROR: no *.v6_clinical.tsv under ${V6_BASE}" >&2
    exit 1
fi

mkdir -p "${REPORT_DIR}"
MANIFEST="${REPORT_DIR}/MANIFEST.tsv"
printf "sample\tichor_pdfs\tichor_params\tv6_clinical\tv6_filtered\tmm_annotated\n" > "${MANIFEST}"

for clin in "${CLIN[@]}"; do
    base="$(basename "${clin}")"
    sample="${base%%.somatic_candidates.withAD.v6_clinical.tsv}"
    # fallback if naming differs (e.g. no .withAD)
    [[ "${sample}" == "${base}" ]] && sample="${base%%.*}"

    echo "[${sample}] collecting"
    s_ichor="${REPORT_DIR}/${sample}/ichorcna"
    s_som="${REPORT_DIR}/${sample}/somatic"
    s_sv="${REPORT_DIR}/${sample}/sv"
    mkdir -p "${s_ichor}" "${s_som}" "${s_sv}"

    # --- ichorCNA genome-wide PDFs + params ---
    n_pdf=0; n_par=0
    if [[ -d "${ICHOR_BASE}/${sample}" ]]; then
        while IFS= read -r f; do cp -f "$f" "${s_ichor}/"; n_pdf=$((n_pdf+1)); done \
            < <(find "${ICHOR_BASE}/${sample}" -name '*_genomeWide*.pdf' 2>/dev/null)
        while IFS= read -r f; do cp -f "$f" "${s_ichor}/"; n_par=$((n_par+1)); done \
            < <(find "${ICHOR_BASE}/${sample}" -name '*.params.txt' 2>/dev/null)
    fi

    # --- v6 somatic reports ---
    v6c="NO"; v6f="NO"
    if [[ -s "${clin}" ]]; then cp -f "${clin}" "${s_som}/"; v6c="YES"; fi
    filt="${clin%.v6_clinical.tsv}.v6_filtered.tsv"
    if [[ -s "${filt}" ]]; then cp -f "${filt}" "${s_som}/"; v6f="YES"; fi

    # --- T2T SV summary ---
    mm="${MM_BASE}/${sample}.mm_annotated.tsv"
    mma="NO"
    if [[ -s "${mm}" ]]; then cp -f "${mm}" "${s_sv}/"; mma="YES"; fi

    printf "%s\t%d\t%d\t%s\t%s\t%s\n" "${sample}" "${n_pdf}" "${n_par}" "${v6c}" "${v6f}" "${mma}" >> "${MANIFEST}"
    echo "    ichor_pdfs=${n_pdf} params=${n_par} v6_clinical=${v6c} v6_filtered=${v6f} mm_annotated=${mma}"
done

echo
echo "Report assembled at: ${REPORT_DIR}"
echo "Manifest:"
cat "${MANIFEST}"
