#!/usr/bin/env bash
#
# STATUS: parked. Wakhan does not produce usable haplotype-specific copy number
# on this data, and the reason is structural rather than a matter of settings.
#
# On the validation cohort the two haplotype coverage tracks were not separated
# anywhere on chr17: the HP1 to HP2 median ratio was 0.76 to 0.87 across 17p,
# across the TP53 panel window, and across 17q alike, where a hemizygous
# deletion at 65% clonality should give roughly 0.2 to 0.3. Regions the tool
# reported as haplotype loss had zero coverage in both haplotypes.
#
# The cause is that phase blocks cannot be chained across the coverage deserts
# between panel windows, so block orientation is effectively arbitrary and
# haplotype-specific depth averages toward the mean. The global purity fit was
# nonetheless accurate (0.6908 against an ichorCNA estimate of 0.693 and FISH
# clone fractions of 65 to 72 per cent), which is why the failure is easy to
# miss: the headline number looks right while the per-segment calls are not.
#
# Uniform coverage is the prerequisite. A shallow whole-genome lane alongside
# the enriched run would supply it; a different tool will not.
#
# For panel-window allelic imbalance use bin/baf_loh_screen.py instead, which
# analyses each window in isolation and does not depend on phase continuity
# between them.
#
#
# run_wakhan_t2t.sh
#
# Allele-specific copy number and CN-LOH detection on T2T using Wakhan.
#
# Rationale
# ---------
# Coverage-based copy number callers cannot detect copy-neutral LOH, and on
# adaptive sampling data the off-target background violates their assumption of
# uniform genome-wide coverage. Wakhan instead segments on phased heterozygous
# SNP allele ratios, which carry allelic information even where read depth is
# thin, and reports allele-specific copy number directly rather than requiring a
# separate BAF track to be joined onto a copy number call.
#
# Wakhan runs in tumour-only mode here, taking a phased Clair3 VCF via
# --tumor-phased-vcf. No whatshap haplotagging stage is required.
#
# The pipeline runs on T2T-CHM13v2.0 rather than hg38. T2T carries no alt or
# decoy contigs, so IGH-side and other repetitive reads are not scattered as
# MAPQ-0 multi-mappers; phase block quality, on which Wakhan's segmentation
# depends, is correspondingly better.
#
# Stages
# ------
#   1. Clair3 genome-wide with phasing on the T2T BAM.
#   2. Wakhan cna, using the phased VCF and the Severus breakpoint VCF.
#
# Stage 1 is the expensive step and has no existing output to reuse: the
# pipeline's own Clair3 phasing runs against hg38.
#
# Usage
# -----
#   ./run_wakhan_t2t.sh <RUN_DIR> <SAMPLE_ID> [MODEL_NAME]
#
# Example
#   ./run_wakhan_t2t.sh results_v7_20260713_24h SAMPLE_ID r1041_e82_400bps_sup_v520
#
# Run detached so the job survives disconnection:
#   setsid nohup ./run_wakhan_t2t.sh <RUN> <SAMPLE> < /dev/null > wakhan.log 2>&1 &

set -o errexit
set -o nounset
set -o pipefail

# ---------------------------------------------------------------------------
# Arguments and paths
# ---------------------------------------------------------------------------

if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <RUN_DIR> <SAMPLE_ID> [MODEL_NAME]" >&2
    exit 1
fi

RUN_DIR="$1"
SAMPLE="$2"
MODEL_NAME="${3:-r1041_e82_400bps_sup_v520}"

REFERENCE=/goast/nikhil_awgs_testing/t2t/refs/chm13v2.0.ucsc.fa
MODEL_DIR="/goast/hemat_data/references/clair3_models/${MODEL_NAME}"
CENTROMERE_BED=/goast/hemat_data/references/T2T/chm13v2.0_censat_v2.1.bed

BAM="${RUN_DIR}/t2t/bams/${SAMPLE}.t2t.bam"
SEVERUS_VCF="${RUN_DIR}/t2t/calls/severus/${SAMPLE}/severus_out/${SAMPLE}.severus.vcf"

CLAIR3_OUT="${RUN_DIR}/t2t/clair3_phased/${SAMPLE}"
WAKHAN_OUT="${RUN_DIR}/t2t/wakhan/${SAMPLE}"

THREADS_CLAIR3=64
THREADS_WAKHAN=16

# Conda is activated inside the script because a detached job does not inherit
# the interactive shell environment.
CONDA_PROFILE=/home/hemat/anaconda3/etc/profile.d/conda.sh

# ---------------------------------------------------------------------------
# Preflight
#
# Every input is checked before any compute is started. Clair3 genome-wide on
# an ONT whole-genome BAM runs for hours, and discovering a missing reference or
# a contig-naming mismatch after that is an expensive way to learn it.
# ---------------------------------------------------------------------------

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

fail() {
    printf '[%s] ERROR: %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2
    exit 1
}

log "Sample ${SAMPLE}, run directory ${RUN_DIR}"

for path in "$BAM" "$REFERENCE" "$CENTROMERE_BED"; do
    [[ -f "$path" ]] || fail "Required input not found: ${path}"
done

[[ -d "$MODEL_DIR" ]] || fail "Clair3 model directory not found: ${MODEL_DIR}"

if [[ ! -f "$SEVERUS_VCF" ]]; then
    log "WARNING: Severus VCF not found at ${SEVERUS_VCF}"
    log "         Wakhan will segment without breakpoint guidance."
    SEVERUS_VCF=""
fi

# ---------------------------------------------------------------------------
# Stage 0: BAM index freshness
#
# A stale index produces silently wrong region queries rather than an error, so
# it is rebuilt whenever it is older than the BAM.
# ---------------------------------------------------------------------------

source "$CONDA_PROFILE"
conda activate awgs_sv

if [[ ! -f "${BAM}.bai" || "${BAM}" -nt "${BAM}.bai" ]]; then
    log "BAM index missing or stale, rebuilding"
    samtools index -@ 16 "$BAM"
fi

# ---------------------------------------------------------------------------
# Contig naming consistency
#
# Wakhan requires the BAM, VCF, reference and --contigs argument to agree on
# chr versus NC_ naming. A mismatch yields empty output rather than an error.
# ---------------------------------------------------------------------------

BAM_CONTIG="$(samtools idxstats "$BAM" | head -1 | cut -f1)"
REF_CONTIG="$(head -1 "${REFERENCE}.fai" | cut -f1)"

log "BAM first contig: ${BAM_CONTIG}"
log "Reference first contig: ${REF_CONTIG}"

if [[ "$BAM_CONTIG" != chr* ]]; then
    fail "BAM contigs are not chr-named (${BAM_CONTIG}). This script assumes chr naming."
fi
if [[ "$REF_CONTIG" != chr* ]]; then
    fail "Reference contigs are not chr-named (${REF_CONTIG}). Use the ucsc-named T2T FASTA."
fi

if [[ -n "$SEVERUS_VCF" ]]; then
    SEVERUS_CONTIG="$(grep -v '^#' "$SEVERUS_VCF" | head -1 | cut -f1 || true)"
    log "Severus first contig: ${SEVERUS_CONTIG:-none}"
    if [[ -n "${SEVERUS_CONTIG}" && "$SEVERUS_CONTIG" != chr* ]]; then
        log "WARNING: Severus VCF is not chr-named; dropping breakpoint input"
        SEVERUS_VCF=""
    fi
fi

# ---------------------------------------------------------------------------
# Stage 1: Clair3 genome-wide with phasing
#
# Phasing is required: Wakhan's segmentation operates on phase blocks, and an
# unphased VCF reduces it to unphased coverage analysis. longphase is used for
# phasing as it is substantially faster than whatshap on ONT data at this scale.
# ---------------------------------------------------------------------------

PHASED_VCF="${CLAIR3_OUT}/phased_merge_output.vcf.gz"

if [[ -f "$PHASED_VCF" ]]; then
    log "Phased VCF already present, skipping Clair3: ${PHASED_VCF}"
else
    log "Stage 1: Clair3 genome-wide phasing (this runs for several hours)"
    mkdir -p "$CLAIR3_OUT"

    # Clair3 runs inside the same container the Nextflow pipeline uses.
    # /goast is bind-mounted at the same path so host paths are valid inside
    # the container and no translation is needed. --user preserves host
    # ownership of the output.
    ABS_BAM="$(readlink -f "$BAM")"
    ABS_REF="$(readlink -f "$REFERENCE")"
    ABS_MODEL="$(readlink -f "$MODEL_DIR")"
    ABS_OUT="$(readlink -f "$CLAIR3_OUT")"

    docker run --rm \
        --user "$(id -u):$(id -g)" \
        -v /goast:/goast \
        -v "${ABS_OUT}:/clair3_out" \
        -w /clair3_out \
        hkubal/clair3:latest \
        /opt/bin/run_clair3.sh \
            --bam_fn="$ABS_BAM" \
            --ref_fn="$ABS_REF" \
            --threads="$THREADS_CLAIR3" \
            --platform=ont \
            --model_path="$ABS_MODEL" \
            --output=/clair3_out \
            --enable_phasing \
            --longphase_for_phasing \
            --include_all_ctgs \
            --sample_name="$SAMPLE"

    for vcf in "${CLAIR3_OUT}/merge_output.vcf.gz" "${CLAIR3_OUT}/phased_merge_output.vcf.gz"; do
        if [[ -s "$vcf" && ! -f "${vcf}.tbi" ]]; then
            tabix -p vcf "$vcf"
        fi
    done

    [[ -f "$PHASED_VCF" ]] || fail "Clair3 completed but no phased VCF at ${PHASED_VCF}"
    log "Stage 1 complete"
fi

# Report phasing quality, since Wakhan's output is only as good as the phase
# blocks it is given.
PHASED_COUNT="$(bcftools view -H -g het "$PHASED_VCF" 2>/dev/null \
    | awk -F'\t' '{n=split($9,k,":"); split($10,v,":");
        for(i=1;i<=n;i++) f[k[i]]=v[i];
        if (index(f["GT"], "|") > 0) c++} END{print c+0}')"
TOTAL_HET="$(bcftools view -H -g het "$PHASED_VCF" 2>/dev/null | wc -l)"
log "Heterozygous sites: ${TOTAL_HET}, of which phased: ${PHASED_COUNT}"

# ---------------------------------------------------------------------------
# Stage 2: Wakhan allele-specific copy number
#
# The purity range is left at its default rather than being constrained to a
# prior from ichorCNA or FISH. An unconstrained fit that independently lands
# near the known tumour fraction is evidence the segmentation is sound; seeding
# it with the expected answer would remove that check.
# ---------------------------------------------------------------------------

log "Stage 2: Wakhan copy number analysis"
conda activate wakhan_env
mkdir -p "$WAKHAN_OUT"

WAKHAN_ARGS=(
    all
    --target-bam "$BAM"
    --reference "$REFERENCE"
    --tumor-phased-vcf "$PHASED_VCF"
    --centromere-bed "$CENTROMERE_BED"
    --out-dir-plots "$WAKHAN_OUT"
    --threads "$THREADS_WAKHAN"
    --contigs "chr1-22,chrX"
    --genome-name "$SAMPLE"
    --reference-name "T2T-CHM13v2.0"
    --pdf-enable
)

if [[ -n "$SEVERUS_VCF" ]]; then
    WAKHAN_ARGS+=(--breakpoints "$SEVERUS_VCF")
    log "Using Severus breakpoints for segmentation"
fi

log "Command: wakhan ${WAKHAN_ARGS[*]}"
wakhan "${WAKHAN_ARGS[@]}"

log "Stage 2 complete. Output in ${WAKHAN_OUT}"
ls -la "$WAKHAN_OUT" || true

log "Done for ${SAMPLE}"
