#!/usr/bin/env bash
#
# compare_t2t_outputs.sh
# ======================
#
# Compare Nextflow T2T-track outputs against production-bash T2T outputs
# for a single sample. Reports:
#
#   * Per-caller (Sniffles, CuteSV, Severus): total / PASS / SVTYPE distribution
#   * Per-caller: BND counts touching the IGH locus (translocation-relevant)
#   * SURVIVOR-merged VCF: total / SUPP_VEC histogram
#   * Merged VCF: BND counts touching IGH + each known partner gene window
#
# This script is read-only. It never modifies the production trees.
#
# Usage:
#   bin/compare_t2t_outputs.sh <sample_id> <bash_root> <nf_root>
#
# Example after the first single-sample Nextflow run:
#   bin/compare_t2t_outputs.sh 11F202612108 \
#       /goast/nikhil_awgs_testing \
#       /goast/nikhil_awgs_testing/mm-awgs-nextflow/results_t2t_only
#
# As a sanity check the script can be pointed at the production tree on both
# sides (bash_root == nf_root); every diff column should then be zero.
#

set -euo pipefail

source /home/hemat/anaconda3/etc/profile.d/conda.sh
conda activate awgs_sv

if [ $# -lt 3 ]; then
    cat >&2 <<EOF
Usage: $0 <sample_id> <bash_root> <nf_root>

  sample_id   E.g. 11F202612108
  bash_root   Root containing t2t/calls/{sniffles,cutesv,severus,merged}/
              For production, this is /goast/nikhil_awgs_testing
  nf_root     Root containing t2t/calls/{sniffles,cutesv,severus,merged}/
              For the Nextflow port, this is .../mm-awgs-nextflow/<outdir>
EOF
    exit 1
fi

SAMPLE="$1"
BASH_ROOT="$2"
NF_ROOT="$3"

# -----------------------------------------------------------------------------
# Region anchors (T2T-CHM13v2.0, chr-named). Coordinates are picked liberally
# enough to bracket the v6 panel windows for these loci so we catch BNDs that
# fell inside the panel but are not exactly at the gene body.
# -----------------------------------------------------------------------------
declare -A REGION_CHR REGION_START REGION_END
REGION_CHR[IGH]=chr14;   REGION_START[IGH]=99000000;   REGION_END[IGH]=101200000
REGION_CHR[IGK]=chr2;    REGION_START[IGK]=88000000;   REGION_END[IGK]=91500000
REGION_CHR[IGL]=chr22;   REGION_START[IGL]=22000000;   REGION_END[IGL]=23800000
REGION_CHR[CCND1]=chr11; REGION_START[CCND1]=69100000; REGION_END[CCND1]=70200000   # t(11;14) BCR
REGION_CHR[NSD2]=chr4;   REGION_START[NSD2]=1500000;   REGION_END[NSD2]=2300000     # t(4;14) BCR
REGION_CHR[MAF]=chr16;   REGION_START[MAF]=78000000;   REGION_END[MAF]=87000000     # WWOX intron 8 to MAF
REGION_CHR[MAFB]=chr20;  REGION_START[MAFB]=41900000;  REGION_END[MAFB]=43000000
REGION_CHR[MYC]=chr8;    REGION_START[MYC]=126000000;  REGION_END[MYC]=131500000

REGIONS_ORDER=(IGH IGK IGL CCND1 NSD2 MAF MAFB MYC)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

# Resolve the Severus VCF inside a given root. Severus' output filename
# varies across versions and somatic/all modes; the Nextflow port also
# publishes under a `severus_out` subdir instead of the per-sample dir
# the production bash uses. Try both layouts.
resolve_severus_vcf() {
    local root="$1"
    local sample="$2"
    # Layout A: bash production       -> <root>/t2t/calls/severus/<sample>/...
    # Layout B: Nextflow v0.1         -> <root>/t2t/calls/severus/severus_out/...
    # Layout C: Nextflow v0.2 cohort  -> <root>/t2t/calls/severus/<sample>/severus_out/...
    local bases=(
        "${root}/t2t/calls/severus/${sample}"
        "${root}/t2t/calls/severus/severus_out"
        "${root}/t2t/calls/severus/${sample}/severus_out"
    )
    for base in "${bases[@]}"; do
        for cand in \
            "${base}/${sample}.severus.vcf"           \
            "${base}/severus_somatic.vcf"             \
            "${base}/severus_all.vcf"                 \
            "${base}/somatic_SVs/severus_somatic.vcf" \
            "${base}/all_SVs/severus_all.vcf"
        do
            if [ -s "$cand" ]; then
                echo "$cand"
                return 0
            fi
        done
    done
    return 1
}

# Count total non-header records in a VCF (handles plain and bgzipped).
count_total() {
    local vcf="$1"
    if [[ "$vcf" == *.gz ]]; then
        bcftools view -H "$vcf" 2>/dev/null | wc -l
    else
        bcftools view -H "$vcf" 2>/dev/null | wc -l
    fi
}

# Count PASS records.
count_pass() {
    local vcf="$1"
    bcftools view -H -f PASS "$vcf" 2>/dev/null | wc -l
}

# Count records with a given SVTYPE.
count_svtype() {
    local vcf="$1"
    local svtype="$2"
    bcftools view -H -i "INFO/SVTYPE=\"${svtype}\"" "$vcf" 2>/dev/null | wc -l
}

# Count BND records where EITHER breakpoint side falls inside the named region.
# For SURVIVOR-merged + Sniffles + CuteSV all the BND ALT format is parsable;
# Severus uses the same VCF 4.2 BND grammar.
count_bnd_in_region() {
    local vcf="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    # The primary-side breakpoint is just CHROM/POS. The mate side is encoded
    # in ALT (e.g. N]chr11:69400000]). We do not need to parse ALT precisely:
    # by symmetry, every cross-chromosome BND appears in the VCF twice (one
    # entry per mate). Counting "CHROM=chr AND start <= POS < end" catches the
    # near-side mate; mate-side appears in the partner region's row. For a
    # single-region count this is exact.
    bcftools view -H -i \
        "INFO/SVTYPE=\"BND\" && CHROM=\"${chr}\" && POS>=${start} && POS<${end}" \
        "$vcf" 2>/dev/null | wc -l
}

# Count BND records spanning a specific (region_a, region_b) pair, i.e. one
# side in region_a AND the other side in region_b. Reads ALT to identify the
# mate side. We count unordered pairs (a partner of either side hits).
count_bnd_spanning_pair() {
    local vcf="$1"
    local chr_a="$2"; local start_a="$3"; local end_a="$4"
    local chr_b="$5"; local start_b="$6"; local end_b="$7"
    # Uses POSIX-portable awk (no match()-with-array, no gensub) so it works
    # under both gawk and mawk. The BND ALT formats supported:
    #   N]chrA:posA]   [chrA:posA[N   N[chrA:posA[   ]chrA:posA]N
    bcftools view -H -i "INFO/SVTYPE=\"BND\"" "$vcf" 2>/dev/null \
        | awk -v ca="$chr_a" -v sa="$start_a" -v ea="$end_a" \
              -v cb="$chr_b" -v sb="$start_b" -v eb="$end_b" '
        function in_region(chrom, pos, c, s, e) {
            return (chrom == c && pos >= s && pos < e)
        }
        {
            chrom = $1
            pos   = $2 + 0
            alt   = $5

            # Find first and second bracket (either [ or ])
            first = 0; second = 0
            n = length(alt)
            for (i = 1; i <= n; i++) {
                ch = substr(alt, i, 1)
                if (ch == "[" || ch == "]") {
                    if (first == 0) first = i
                    else { second = i; break }
                }
            }
            if (first == 0 || second == 0) next

            inner = substr(alt, first + 1, second - first - 1)
            np = split(inner, parts, ":")
            if (np < 2) next
            mate_chrom = parts[1]
            mate_pos   = parts[2] + 0

            # Match either ordering of the (a, b) pair
            if (in_region(chrom, pos, ca, sa, ea) && in_region(mate_chrom, mate_pos, cb, sb, eb)) count++
            else if (in_region(chrom, pos, cb, sb, eb) && in_region(mate_chrom, mate_pos, ca, sa, ea)) count++
        }
        END { print (count ? count : 0) }
        '
}

# SUPP_VEC histogram. Output: lines of "<bits> <count>" sorted by bits.
supp_vec_histogram() {
    local vcf="$1"
    bcftools view -H "$vcf" 2>/dev/null \
        | grep -oP 'SUPP_VEC=\K[01]+' \
        | sort \
        | uniq -c \
        | awk '{print $2, $1}'
}

# Print a row of a comparison table. cols are: label bash_value nf_value
print_row() {
    local label="$1"
    local bash_v="$2"
    local nf_v="$3"
    local diff=$(( nf_v - bash_v ))
    local marker=""
    if [ "$diff" -ne 0 ]; then
        marker="  <-- diff=${diff}"
    fi
    printf "  %-32s  %12s  %12s  %s\n" "$label" "$bash_v" "$nf_v" "$marker"
}

# Banner / header for a comparison block.
banner() {
    echo ""
    echo "================================================================"
    echo "$*"
    echo "================================================================"
}

# -----------------------------------------------------------------------------
# Locate VCFs in each tree
# -----------------------------------------------------------------------------
echo "================================================================"
echo "Comparing T2T outputs for sample: ${SAMPLE}"
echo "  bash root:     ${BASH_ROOT}"
echo "  nextflow root: ${NF_ROOT}"
echo "================================================================"

# Per-caller paths
BASH_SNIFFLES="${BASH_ROOT}/t2t/calls/sniffles/${SAMPLE}.sniffles.t2t.vcf.gz"
NF_SNIFFLES="${NF_ROOT}/t2t/calls/sniffles/${SAMPLE}.sniffles.t2t.vcf.gz"
BASH_CUTESV="${BASH_ROOT}/t2t/calls/cutesv/${SAMPLE}.cutesv.t2t.vcf.gz"
NF_CUTESV="${NF_ROOT}/t2t/calls/cutesv/${SAMPLE}.cutesv.t2t.vcf.gz"
BASH_MERGED="${BASH_ROOT}/t2t/calls/merged/${SAMPLE}.merged.vcf.gz"
NF_MERGED="${NF_ROOT}/t2t/calls/merged/${SAMPLE}.merged.vcf.gz"

BASH_SEVERUS=$(resolve_severus_vcf "$BASH_ROOT" "$SAMPLE" || true)
NF_SEVERUS=$(resolve_severus_vcf "$NF_ROOT" "$SAMPLE" || true)

# Existence guard
missing=0
for f in "$BASH_SNIFFLES" "$NF_SNIFFLES" "$BASH_CUTESV" "$NF_CUTESV" \
         "$BASH_MERGED" "$NF_MERGED"
do
    if [ ! -s "$f" ]; then
        echo "MISSING:  $f" >&2
        missing=$((missing + 1))
    fi
done
if [ -z "$BASH_SEVERUS" ]; then
    echo "MISSING:  Severus VCF for ${SAMPLE} under ${BASH_ROOT}/t2t/calls/severus/${SAMPLE}/" >&2
    missing=$((missing + 1))
fi
if [ -z "$NF_SEVERUS" ]; then
    echo "MISSING:  Severus VCF for ${SAMPLE} under ${NF_ROOT}/t2t/calls/severus/${SAMPLE}/" >&2
    missing=$((missing + 1))
fi
if [ "$missing" -gt 0 ]; then
    echo "" >&2
    echo "Cannot continue; $missing input files not found." >&2
    exit 2
fi

# -----------------------------------------------------------------------------
# Per-caller comparison helper
# -----------------------------------------------------------------------------
compare_caller() {
    local caller="$1"
    local bash_vcf="$2"
    local nf_vcf="$3"

    banner "Caller: ${caller}"
    printf "  %-32s  %12s  %12s\n" "" "BASH" "NEXTFLOW"

    local bash_total=$(count_total "$bash_vcf")
    local nf_total=$(  count_total "$nf_vcf")
    print_row "total records" "$bash_total" "$nf_total"

    local bash_pass=$(count_pass "$bash_vcf")
    local nf_pass=$(  count_pass "$nf_vcf")
    print_row "PASS records"  "$bash_pass" "$nf_pass"

    for sv in BND DEL DUP INV INS; do
        local bash_n=$(count_svtype "$bash_vcf" "$sv")
        local nf_n=$(  count_svtype "$nf_vcf"   "$sv")
        print_row "  SVTYPE=${sv}" "$bash_n" "$nf_n"
    done

    echo ""
    echo "  BND counts with near-side breakpoint in each region:"
    printf "  %-32s  %12s  %12s\n" "" "BASH" "NEXTFLOW"
    for r in "${REGIONS_ORDER[@]}"; do
        local bash_n=$(count_bnd_in_region "$bash_vcf" \
                        "${REGION_CHR[$r]}" "${REGION_START[$r]}" "${REGION_END[$r]}")
        local nf_n=$(  count_bnd_in_region "$nf_vcf"   \
                        "${REGION_CHR[$r]}" "${REGION_START[$r]}" "${REGION_END[$r]}")
        print_row "  $r (${REGION_CHR[$r]}:${REGION_START[$r]}-${REGION_END[$r]})" \
                  "$bash_n" "$nf_n"
    done

    echo ""
    echo "  BNDs spanning canonical-partner pairs (both sides hit):"
    printf "  %-32s  %12s  %12s\n" "" "BASH" "NEXTFLOW"
    for partner in CCND1 NSD2 MAF MAFB MYC; do
        local bash_n=$(count_bnd_spanning_pair "$bash_vcf" \
                        "${REGION_CHR[IGH]}" "${REGION_START[IGH]}" "${REGION_END[IGH]}" \
                        "${REGION_CHR[$partner]}" "${REGION_START[$partner]}" "${REGION_END[$partner]}")
        local nf_n=$(  count_bnd_spanning_pair "$nf_vcf" \
                        "${REGION_CHR[IGH]}" "${REGION_START[IGH]}" "${REGION_END[IGH]}" \
                        "${REGION_CHR[$partner]}" "${REGION_START[$partner]}" "${REGION_END[$partner]}")
        print_row "  IGH<->$partner" "$bash_n" "$nf_n"
    done
}

# -----------------------------------------------------------------------------
# Run per-caller comparisons
# -----------------------------------------------------------------------------
compare_caller "Sniffles2"   "$BASH_SNIFFLES" "$NF_SNIFFLES"
compare_caller "CuteSV"      "$BASH_CUTESV"   "$NF_CUTESV"
compare_caller "Severus"     "$BASH_SEVERUS"  "$NF_SEVERUS"

# -----------------------------------------------------------------------------
# Merged VCF: total + SUPP_VEC histogram + canonical partner BNDs
# -----------------------------------------------------------------------------
banner "Merged VCF (SURVIVOR ensemble)"
printf "  %-32s  %12s  %12s\n" "" "BASH" "NEXTFLOW"

bash_total=$(count_total "$BASH_MERGED")
nf_total=$(  count_total "$NF_MERGED")
print_row "total merged records" "$bash_total" "$nf_total"

echo ""
echo "  SUPP_VEC histogram (bit order: Sniffles, CuteSV, Severus):"
printf "  %-32s  %12s  %12s\n" "" "BASH" "NEXTFLOW"

# Collect histograms into associative arrays
declare -A bash_hist nf_hist
while read -r bits cnt; do bash_hist[$bits]=$cnt; done < <(supp_vec_histogram "$BASH_MERGED")
while read -r bits cnt; do   nf_hist[$bits]=$cnt; done < <(supp_vec_histogram "$NF_MERGED")

# Union of bit patterns, sorted for stable output
all_bits=$( { printf '%s\n' "${!bash_hist[@]}"; printf '%s\n' "${!nf_hist[@]}"; } \
             | sort -u )
for b in $all_bits; do
    bv=${bash_hist[$b]:-0}
    nv=${nf_hist[$b]:-0}
    case "$b" in
        100) lbl="Sniffles only" ;;
        010) lbl="CuteSV only" ;;
        001) lbl="Severus only" ;;
        110) lbl="Sniffles + CuteSV" ;;
        101) lbl="Sniffles + Severus" ;;
        011) lbl="CuteSV + Severus" ;;
        111) lbl="all three (3-caller agreement)" ;;
        *)   lbl="$b (unexpected)" ;;
    esac
    print_row "  $b  $lbl" "$bv" "$nv"
done

echo ""
echo "  Merged-VCF BNDs spanning canonical pairs:"
printf "  %-32s  %12s  %12s\n" "" "BASH" "NEXTFLOW"
for partner in CCND1 NSD2 MAF MAFB MYC; do
    bash_n=$(count_bnd_spanning_pair "$BASH_MERGED" \
              "${REGION_CHR[IGH]}" "${REGION_START[IGH]}" "${REGION_END[IGH]}" \
              "${REGION_CHR[$partner]}" "${REGION_START[$partner]}" "${REGION_END[$partner]}")
    nf_n=$(  count_bnd_spanning_pair "$NF_MERGED" \
              "${REGION_CHR[IGH]}" "${REGION_START[IGH]}" "${REGION_END[IGH]}" \
              "${REGION_CHR[$partner]}" "${REGION_START[$partner]}" "${REGION_END[$partner]}")
    print_row "  IGH<->$partner" "$bash_n" "$nf_n"
done

echo ""
echo "================================================================"
echo "Done."
echo "================================================================"
echo ""
echo "How to read this:"
echo "  - 'BASH' column = production-bash output that already lives on disk."
echo "  - 'NEXTFLOW' column = current Nextflow port output."
echo "  - 'diff=N' annotation flags any count where the two disagree."
echo ""
echo "Acceptance criteria for the v0.1.1 single-sample validation:"
echo "  - PASS counts per caller match exactly."
echo "  - SUPP_VEC histogram pattern proportions match (small absolute"
echo "    differences acceptable if the BAM was re-realigned post-fix)."
echo "  - IGH-spanning BND counts match the FISH-known translocation for"
echo "    this sample (sample 11F202612108 has t(11;14) 95% by FISH, so"
echo "    IGH<->CCND1 should be > 0 in both columns)."
