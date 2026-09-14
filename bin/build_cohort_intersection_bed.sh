#!/usr/bin/env bash
#
# build_cohort_intersection_bed.sh
#
# Build the intersection of the v4, v5, v6 and v7 capture panels, for use by
# cohort-level processes that compare every sample against a shared baseline.
#
# Why this is needed
# ------------------
# BAF_LOH_SCREEN and BAF_CN_PLOTS compare each region against the cohort median
# for that same region. They take one BED for all samples, which is correct for
# a cohort statistic and wrong if that BED is any single panel version: the
# cohort spans v4, v5, v6 and v7, so a v7 BED asks six samples to contribute
# off-target depth at regions the other fifteen captured on-target. That biases
# the median at exactly the regions v7 added, which is where the difference
# between panel versions lives.
#
# In the intersection every sample has the same capture status at every region,
# so the median means one thing. The cost is coverage: regions present in only
# some panels are dropped from the cohort screen. They are still analysed per
# sample, against that sample's own panel.
#
# This is deliberately separate from params.panel_bed_hg38, which remains the
# per-sample default. The two uses have diverged and should not share a
# parameter.
#
# Usage:
#   bash build_cohort_intersection_bed.sh [assets_dir]
#
# Writes assets/aWGS_cohort_intersection_hg38.bed and prints the size of each
# input panel alongside the intersection, so the cost is visible before the
# result is used.

set -euo pipefail

source /home/hemat/anaconda3/etc/profile.d/conda.sh
conda activate awgs_sv

ASSETS="${1:-/goast/mm-awgs-nextflow/assets}"
OUT="${ASSETS}/aWGS_cohort_intersection_hg38.bed"

PANELS=(
    "${ASSETS}/aWGS_MMfocused_v4_hg38.bed"
    "${ASSETS}/aWGS_MMfocused_v5_hg38.bed"
    "${ASSETS}/aWGS_MMfocused_v6_hg38.bed"
    "${ASSETS}/aWGS_PCN_v7_hg38.bed"
)

for bed in "${PANELS[@]}"; do
    if [[ ! -f "$bed" ]]; then
        echo "ERROR: missing panel BED: $bed" >&2
        exit 1
    fi
done

total_bp() {
    awk 'BEGIN{s=0} !/^(#|track|browser)/ && NF>=3 {s += $3 - $2} END{printf "%d", s}' "$1"
}

region_count() {
    awk 'BEGIN{n=0} !/^(#|track|browser)/ && NF>=3 {n++} END{printf "%d", n}' "$1"
}

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

echo "Input panels"
for bed in "${PANELS[@]}"; do
    printf "  %-40s %3d regions  %12s bp\n" \
        "$(basename "$bed")" "$(region_count "$bed")" \
        "$(printf "%'d" "$(total_bp "$bed")")"
done

# Sort and merge each panel first. Intersecting unsorted or self-overlapping
# intervals gives duplicate output rows and an inflated total.
for i in "${!PANELS[@]}"; do
    sort -k1,1 -k2,2n "${PANELS[$i]}" \
        | bedtools merge -i - > "${TMP}/panel_${i}.bed"
done

cp "${TMP}/panel_0.bed" "${TMP}/acc.bed"
for i in 1 2 3; do
    bedtools intersect -a "${TMP}/acc.bed" -b "${TMP}/panel_${i}.bed" \
        > "${TMP}/next.bed"
    mv "${TMP}/next.bed" "${TMP}/acc.bed"
done

sort -k1,1 -k2,2n "${TMP}/acc.bed" | bedtools merge -i - > "$OUT"

echo ""
echo "Intersection"
printf "  %-40s %3d regions  %12s bp\n" \
    "$(basename "$OUT")" "$(region_count "$OUT")" \
    "$(printf "%'d" "$(total_bp "$OUT")")"

v7_bp="$(total_bp "${PANELS[3]}")"
int_bp="$(total_bp "$OUT")"
if [[ "$v7_bp" -gt 0 ]]; then
    echo ""
    awk -v a="$int_bp" -v b="$v7_bp" \
        'BEGIN{printf "  retains %.1f%% of the v7 panel\n", 100*a/b}'
fi

echo ""
echo "Written: $OUT"
echo "Point params.cohort_bed_hg38 at it."
echo ""
echo "Note: the intersection has no gene names. The cohort screen uses it for"
echo "region membership only, not for annotation, so that is fine; do not"
echo "reuse this file anywhere that expects a name column."
