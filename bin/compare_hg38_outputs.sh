#!/usr/bin/env bash
#
# compare_hg38_outputs.sh
#
# Validate the Nextflow hg38 track against the production-bash hg38 outputs
# for a single sample. Mirrors compare_t2t_outputs.sh in spirit.
#
# Modules and how they are checked:
#   CLAIR3_PHASED         - compared against production bash (total / PASS /
#                           phased-het counts on merge_output and
#                           phased_merge_output).
#   VEP_ANNOTATE_CLAIR3   - compared against production bash (all-annotated
#                           rows and somatic-candidate rows).
#   CLAIRS_TO             - BLIND. Production bash did not run ClairS-TO for the
#                           18h snapshot sample, so we only print the Nextflow
#                           SNV/indel counts as a sanity check (no comparison).
#   ICHORCNA              - BLIND. Same reason. We print tumor fraction and
#                           ploidy for interpretation only.
#
# Usage:
#   compare_hg38_outputs.sh <sample_id> <bash_root> <nextflow_root>
#
# Example:
#   compare_hg38_outputs.sh 11F202612108_18h \
#       /goast/nikhil_awgs_testing \
#       /goast/mm-awgs-nextflow/results_hg38_only
#
# A line annotated with "<-- diff=N" flags a count where the two pipelines
# disagree. Grep for that marker to count real disagreements:
#   grep -c "<-- diff=" output.txt
#
set -uo pipefail

SAMPLE="${1:?need sample_id}"
BASH_ROOT="${2:?need bash_root}"
NF_ROOT="${3:?need nextflow_root}"

# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------

# Count total records in a VCF (handles .gz and plain).
vcf_total() {
    local vcf="$1"
    [ -s "$vcf" ] || { echo "NA"; return; }
    bcftools view -H "$vcf" 2>/dev/null | wc -l | tr -d ' '
}

# Count PASS records in a VCF.
vcf_pass() {
    local vcf="$1"
    [ -s "$vcf" ] || { echo "NA"; return; }
    bcftools view -f PASS -H "$vcf" 2>/dev/null | wc -l | tr -d ' '
}

# Count phased genotypes (GT containing a pipe) in a VCF.
vcf_phased() {
    local vcf="$1"
    [ -s "$vcf" ] || { echo "NA"; return; }
    # Pull the GT field of the first sample column and count those with a pipe.
    bcftools query -f '[%GT]\n' "$vcf" 2>/dev/null | grep -c '|'
}

# Count data rows (excluding header) in a TSV.
tsv_rows() {
    local tsv="$1"
    [ -s "$tsv" ] || { echo "NA"; return; }
    local n
    n=$(wc -l < "$tsv" | tr -d ' ')
    echo $(( n > 0 ? n - 1 : 0 ))
}

# Print a labelled comparison row, flagging diffs.
# Args: label bash_val nf_val
cmp_row() {
    local label="$1" b="$2" n="$3"
    local flag=""
    if [ "$b" != "NA" ] && [ "$n" != "NA" ] && [ "$b" != "$n" ]; then
        flag="   <-- diff=$(( n - b ))"
    fi
    printf "  %-44s %10s %12s%s\n" "$label" "$b" "$n" "$flag"
}

# Print a single blind metric (Nextflow-only, no bash comparison).
blind_row() {
    local label="$1" n="$2"
    printf "  %-44s %10s %12s\n" "$label" "(blind)" "$n"
}

# Resolve the Clair3 output dir under a root. Production bash uses
# <root>/hg38/clair3_phased/<sample>/; the Nextflow port publishes into
# <root>/hg38/clair3_phased/clair3_out/.
resolve_clair3_dir() {
    local root="$1" sample="$2"
    for d in \
        "${root}/hg38/clair3_phased/${sample}" \
        "${root}/hg38/clair3_phased/clair3_out"
    do
        if [ -s "${d}/merge_output.vcf.gz" ]; then
            echo "$d"
            return 0
        fi
    done
    return 1
}

# Resolve the VEP-annotated dir under a root. Production bash uses
# <root>/hg38/calls/annotated_clair3/<sample>/; the Nextflow port publishes
# into <root>/hg38/calls/annotated_clair3/vep_out/.
resolve_vep_dir() {
    local root="$1" sample="$2"
    for d in \
        "${root}/hg38/calls/annotated_clair3/${sample}" \
        "${root}/hg38/calls/annotated_clair3/vep_out"
    do
        if [ -s "${d}/${sample}.all_annotated.tsv" ] || [ -s "${d}/${sample}.somatic_candidates.tsv" ]; then
            echo "$d"
            return 0
        fi
    done
    return 1
}

# Resolve the ClairS-TO output dir (Nextflow only; bash absent for 18h).
resolve_clairs_to_dir() {
    local root="$1" sample="$2"
    for d in \
        "${root}/hg38/calls/clairs_to/${sample}" \
        "${root}/hg38/calls/clairs_to/clairs_to_out"
    do
        if [ -d "$d" ]; then
            echo "$d"
            return 0
        fi
    done
    return 1
}

# Resolve the ichorCNA params file (Nextflow only; bash absent for 18h).
resolve_ichorcna_params() {
    local root="$1" sample="$2"
    for p in \
        "${root}/hg38/calls/ichorcna/${sample}/${sample}.params.txt" \
        "${root}/hg38/calls/ichorcna/ichorcna_out/${sample}.params.txt" \
        "${root}/hg38/calls/ichorcna/${sample}.params.txt"
    do
        if [ -s "$p" ]; then
            echo "$p"
            return 0
        fi
    done
    return 1
}

echo "================================================================"
echo "Comparing hg38 outputs for sample: ${SAMPLE}"
echo "  bash root:     ${BASH_ROOT}"
echo "  nextflow root: ${NF_ROOT}"
echo "================================================================"

# ----------------------------------------------------------------------------
# CLAIR3_PHASED  (compared against production bash)
# ----------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "CLAIR3_PHASED  (validated against production bash)"
echo "================================================================"

bash_c3=$(resolve_clair3_dir "$BASH_ROOT" "$SAMPLE" || true)
nf_c3=$(resolve_clair3_dir "$NF_ROOT" "$SAMPLE" || true)

if [ -z "$bash_c3" ]; then echo "  WARNING: production-bash Clair3 dir not found"; fi
if [ -z "$nf_c3" ];   then echo "  WARNING: nextflow Clair3 dir not found"; fi

printf "  %-44s %10s %12s\n" "" "BASH" "NEXTFLOW"

# merge_output
b_merge="${bash_c3:+$bash_c3/merge_output.vcf.gz}"
n_merge="${nf_c3:+$nf_c3/merge_output.vcf.gz}"
cmp_row "merge_output total"        "$(vcf_total "${b_merge:-}")"  "$(vcf_total "${n_merge:-}")"
cmp_row "merge_output PASS"         "$(vcf_pass  "${b_merge:-}")"  "$(vcf_pass  "${n_merge:-}")"

# phased_merge_output
b_phase="${bash_c3:+$bash_c3/phased_merge_output.vcf.gz}"
n_phase="${nf_c3:+$nf_c3/phased_merge_output.vcf.gz}"
cmp_row "phased_merge_output total" "$(vcf_total  "${b_phase:-}")" "$(vcf_total  "${n_phase:-}")"
cmp_row "phased_merge_output PASS"  "$(vcf_pass   "${b_phase:-}")" "$(vcf_pass   "${n_phase:-}")"
cmp_row "phased_merge_output phased GTs" "$(vcf_phased "${b_phase:-}")" "$(vcf_phased "${n_phase:-}")"

# ----------------------------------------------------------------------------
# VEP_ANNOTATE_CLAIR3  (compared against production bash)
# ----------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "VEP_ANNOTATE_CLAIR3  (validated against production bash)"
echo "================================================================"

bash_vep=$(resolve_vep_dir "$BASH_ROOT" "$SAMPLE" || true)
nf_vep=$(resolve_vep_dir "$NF_ROOT" "$SAMPLE" || true)

if [ -z "$bash_vep" ]; then echo "  WARNING: production-bash VEP dir not found"; fi
if [ -z "$nf_vep" ];   then echo "  WARNING: nextflow VEP dir not found"; fi

printf "  %-44s %10s %12s\n" "" "BASH" "NEXTFLOW"

b_all="${bash_vep:+$bash_vep/${SAMPLE}.all_annotated.tsv}"
n_all="${nf_vep:+$nf_vep/${SAMPLE}.all_annotated.tsv}"
cmp_row "all_annotated rows"        "$(tsv_rows "${b_all:-}")"  "$(tsv_rows "${n_all:-}")"

b_cand="${bash_vep:+$bash_vep/${SAMPLE}.somatic_candidates.tsv}"
n_cand="${nf_vep:+$nf_vep/${SAMPLE}.somatic_candidates.tsv}"
cmp_row "somatic_candidates rows"   "$(tsv_rows "${b_cand:-}")" "$(tsv_rows "${n_cand:-}")"

# PASS VCF record count, if both produced one
b_pass_vcf="${bash_vep:+$bash_vep/${SAMPLE}.pass.vcf.gz}"
n_pass_vcf="${nf_vep:+$nf_vep/${SAMPLE}.pass.vcf.gz}"
cmp_row "pass.vcf records"          "$(vcf_total "${b_pass_vcf:-}")" "$(vcf_total "${n_pass_vcf:-}")"

# ----------------------------------------------------------------------------
# CLAIRS_TO  (blind; production bash absent for 18h)
# ----------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "CLAIRS_TO  (blind sanity-check; no production-bash reference)"
echo "================================================================"

nf_cs=$(resolve_clairs_to_dir "$NF_ROOT" "$SAMPLE" || true)
printf "  %-44s %10s %12s\n" "" "BASH" "NEXTFLOW"
if [ -n "$nf_cs" ]; then
    snv=$(vcf_total "$nf_cs/snv_${SAMPLE}.vcf.gz")
    indel=$(vcf_total "$nf_cs/indel_${SAMPLE}.vcf.gz")
    blind_row "SNV records"   "$snv"
    blind_row "indel records" "$indel"
else
    echo "  WARNING: nextflow ClairS-TO dir not found"
fi

# ----------------------------------------------------------------------------
# ICHORCNA  (blind; production bash absent for 18h)
# ----------------------------------------------------------------------------
echo ""
echo "================================================================"
echo "ICHORCNA  (blind interpretation; no production-bash reference)"
echo "================================================================"

nf_ich=$(resolve_ichorcna_params "$NF_ROOT" "$SAMPLE" || true)
if [ -n "$nf_ich" ]; then
    tf=$(grep -i "Tumor Fraction:" "$nf_ich" | awk -F: '{gsub(/ /,"",$2); print $2}')
    pl=$(grep -i "Ploidy:"         "$nf_ich" | awk -F: '{gsub(/ /,"",$2); print $2}')
    blind_row "tumor fraction" "${tf:-NA}"
    blind_row "ploidy"         "${pl:-NA}"
    echo ""
    echo "  Note: at 18h adaptive-sampling depth, off-target coverage is sparse;"
    echo "  a TF near 0 is expected and not a pipeline error. Rerun on full-depth"
    echo "  data for a meaningful tumor fraction estimate."
else
    echo "  WARNING: nextflow ichorCNA params.txt not found"
fi

echo ""
echo "================================================================"
echo "Done."
echo "================================================================"
echo ""
echo "How to read this:"
echo "  - CLAIR3_PHASED and VEP rows compare BASH vs NEXTFLOW; a 'diff=N'"
echo "    annotation flags any count that disagrees."
echo "  - CLAIRS_TO and ICHORCNA are blind (production bash did not run them"
echo "    for the 18h snapshot), so only Nextflow values are shown."
echo ""
echo "Acceptance for the hg38 single-sample validation:"
echo "  - Clair3 merge_output total and PASS match production bash."
echo "  - phased_merge_output phased-GT count matches (phasing reproduced)."
echo "  - VEP all_annotated and somatic_candidates row counts match."
echo "  - ClairS-TO produces a plausible SNV count (~10-15k expected)."
echo "  - ichorCNA runs to completion (TF interpretation noted above)."
