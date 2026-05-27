#!/usr/bin/env bash
# qc_v6_cohort.sh
# =================
# Regression QC for v6 panel against cohort BAMs.
#
# What this does:
#   1. Auto-detects all *.t2t.bam files in $T2T_BAM_DIR.
#   2. Classifies each BAM as either "full" (final run) or "18h" (early
#      snapshot) based on filename suffix (_18h).
#   3. Runs mosdepth with the v6 chr BED for each BAM (one per-sample
#      output regardless of group).
#   4. Emits TWO cohort coverage matrices to avoid double-counting the
#      same biological sample at two timepoints:
#        v6_cohort_coverage_full.tsv   (full-run BAMs)
#        v6_cohort_coverage_18h.tsv    (18h snapshot BAMs, if any exist)
#   5. If both groups have at least one shared biological sample,
#      emits a paired 18h-vs-full comparison TSV for monitoring
#      adaptive-sampling efficiency over time.
#
# This script is variant-agnostic and sample-agnostic. The only filename
# convention it relies on is the "_18h" suffix that distinguishes early
# snapshots from full-run BAMs.
#
# Usage:
#   ./qc_v6_cohort.sh                       # uses defaults from paths.sh
#   ./qc_v6_cohort.sh /custom/outdir        # override output directory
#
# Requires conda env awgs_sv (mosdepth, samtools).

set -euo pipefail

# Source conda for nohup-safe execution
source /home/hemat/anaconda3/etc/profile.d/conda.sh
conda activate awgs_sv

# Source project paths
PATHS_SH="/goast/nikhil_awgs_testing/config/paths.sh"
if [[ -f "$PATHS_SH" ]]; then
    # shellcheck disable=SC1090
    source "$PATHS_SH"
else
    echo "ERROR: paths.sh not found at $PATHS_SH" >&2
    exit 1
fi

# Resolve inputs
BAM_DIR="${T2T_BAM_DIR:-/goast/nikhil_awgs_testing/t2t/bams}"
PANEL_BED="${T2T_PANEL_BED:-/goast/nikhil_awgs_testing/t2t/beds/aWGS_MMfocused_v6_t2t_chr.bed}"
OUTDIR="${1:-${T2T_QC_DIR:-/goast/nikhil_awgs_testing/t2t/qc}/v6_cohort}"
THREADS="${THREADS:-4}"

mkdir -p "$OUTDIR"

echo "=== v6 cohort regression QC ==="
echo "BAM directory:  $BAM_DIR"
echo "Panel BED:      $PANEL_BED"
echo "Output:         $OUTDIR"
echo "Threads:        $THREADS"
echo

# Sanity checks
[[ -d "$BAM_DIR" ]]    || { echo "ERROR: BAM_DIR not found: $BAM_DIR"; exit 1; }
[[ -f "$PANEL_BED" ]]  || { echo "ERROR: PANEL_BED not found: $PANEL_BED"; exit 1; }

PANEL_MD5=$(md5sum "$PANEL_BED" | cut -d' ' -f1)
EXPECTED_MD5="b9ce72ca3b1ba1294b8d49bad0b7dab2"
if [[ "$PANEL_MD5" != "$EXPECTED_MD5" ]]; then
    echo "WARNING: panel BED MD5 does not match expected v6 chr MD5."
    echo "  Got:      $PANEL_MD5"
    echo "  Expected: $EXPECTED_MD5"
fi

# Discover and classify cohort BAMs
mapfile -t ALL_BAMS < <(find "$BAM_DIR" -maxdepth 1 -type f -name '*.t2t.bam' | sort)
if [[ ${#ALL_BAMS[@]} -eq 0 ]]; then
    echo "ERROR: no *.t2t.bam files found in $BAM_DIR"
    exit 1
fi

FULL_BAMS=()
H18_BAMS=()
for bam in "${ALL_BAMS[@]}"; do
    name=$(basename "$bam" .t2t.bam)
    if [[ "$name" == *_18h ]]; then
        H18_BAMS+=("$bam")
    else
        FULL_BAMS+=("$bam")
    fi
done

echo "Found ${#ALL_BAMS[@]} BAM(s):"
echo "  Full-run BAMs (n=${#FULL_BAMS[@]}):"
for bam in "${FULL_BAMS[@]}"; do echo "    $(basename "$bam")"; done
echo "  18h snapshot BAMs (n=${#H18_BAMS[@]}):"
for bam in "${H18_BAMS[@]}"; do echo "    $(basename "$bam")"; done
echo

# Run mosdepth on each BAM (regardless of group)
PER_SAMPLE_DIR="$OUTDIR/per_sample"
mkdir -p "$PER_SAMPLE_DIR"

run_mosdepth() {
    local bam="$1"
    local sample
    sample=$(basename "$bam" .t2t.bam)
    local prefix="$PER_SAMPLE_DIR/$sample"
    # Skip if mosdepth output already exists and is newer than the BAM
    if [[ -f "${prefix}.regions.bed.gz" && "${prefix}.regions.bed.gz" -nt "$bam" ]]; then
        echo "  skip (cached): $sample"
        return
    fi
    echo "  mosdepth: $sample"
    mosdepth \
        --no-per-base \
        --by "$PANEL_BED" \
        --mapq 0 \
        --threads "$THREADS" \
        "$prefix" \
        "$bam"
}

echo "Computing per-region depth (mosdepth):"
for bam in "${ALL_BAMS[@]}"; do
    run_mosdepth "$bam"
done
echo

# Build cohort coverage matrices using a Python helper.
# Two outputs: one for full-run, one for 18h snapshots.
build_cohort_matrix() {
    local label="$1"; shift
    local out_tsv="$1"; shift
    local bams=("$@")
    if [[ ${#bams[@]} -eq 0 ]]; then
        echo "  $label: no BAMs in group, skipping matrix"
        return
    fi
    echo "  $label: assembling matrix from ${#bams[@]} BAM(s) -> $out_tsv"
    python3 - "$PER_SAMPLE_DIR" "$out_tsv" "${bams[@]}" <<'PYEOF'
import sys, gzip, statistics
from pathlib import Path

per_sample_dir = Path(sys.argv[1])
out_tsv = Path(sys.argv[2])
bams = sys.argv[3:]
sample_ids = [Path(b).name.replace(".t2t.bam", "") for b in bams]

per_sample = {}
for s in sample_ids:
    regions_gz = per_sample_dir / f"{s}.regions.bed.gz"
    if not regions_gz.exists():
        sys.stderr.write(f"missing: {regions_gz}\n")
        sys.exit(1)
    rows = {}
    with gzip.open(regions_gz, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            chrom, start, end, name, depth = parts[0], parts[1], parts[2], parts[3], parts[4]
            rows[(chrom, start, end, name)] = float(depth)
    per_sample[s] = rows

# Use the first sample as the canonical region order (mosdepth preserves BED order).
regions = list(per_sample[sample_ids[0]].keys())

with open(out_tsv, "w") as out:
    # Header
    out.write("chrom\tstart\tend\tname")
    for s in sample_ids:
        out.write(f"\t{s}")
    out.write("\tmean\tmin\tmax\n")
    for r in regions:
        depths = [per_sample[s].get(r, float("nan")) for s in sample_ids]
        valid = [d for d in depths if d == d]
        mean_c = statistics.mean(valid) if valid else float("nan")
        min_c = min(valid) if valid else float("nan")
        max_c = max(valid) if valid else float("nan")
        out.write("\t".join([
            r[0], r[1], r[2], r[3],
            *[f"{d:.2f}" for d in depths],
            f"{mean_c:.2f}", f"{min_c:.2f}", f"{max_c:.2f}",
        ]) + "\n")
PYEOF
}

echo "Building cohort matrices:"
COHORT_FULL_TSV="$OUTDIR/v6_cohort_coverage_full.tsv"
COHORT_18H_TSV="$OUTDIR/v6_cohort_coverage_18h.tsv"
build_cohort_matrix "full" "$COHORT_FULL_TSV" "${FULL_BAMS[@]}"
build_cohort_matrix "18h"  "$COHORT_18H_TSV"  "${H18_BAMS[@]}"
echo

# Paired comparison: where the same biological sample has both 18h and full,
# emit a TSV of per-region depth at 18h vs full and the ratio (full/18h).
PAIRED_TSV="$OUTDIR/v6_18h_vs_full.tsv"
if [[ ${#H18_BAMS[@]} -gt 0 && ${#FULL_BAMS[@]} -gt 0 ]]; then
    echo "Building paired 18h-vs-full comparison: $PAIRED_TSV"
    python3 - "$PER_SAMPLE_DIR" "$PAIRED_TSV" "${H18_BAMS[@]}" --break "${FULL_BAMS[@]}" <<'PYEOF'
import sys, gzip
from pathlib import Path

argv = sys.argv[1:]
per_sample_dir = Path(argv[0])
out_tsv = Path(argv[1])
rest = argv[2:]
brk = rest.index("--break")
h18_bams = rest[:brk]
full_bams = rest[brk+1:]

def name_of(bam): return Path(bam).name.replace(".t2t.bam", "")

h18_ids = {name_of(b)[:-len("_18h")]: name_of(b) for b in h18_bams}
full_ids = {name_of(b): name_of(b) for b in full_bams}
shared = sorted(set(h18_ids) & set(full_ids))

if not shared:
    sys.stderr.write("No biological samples are present in both groups; skipping paired comparison.\n")
    out_tsv.write_text("")  # empty file as a flag
    sys.exit(0)

def load(s):
    p = per_sample_dir / f"{s}.regions.bed.gz"
    rows = {}
    with gzip.open(p, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 5:
                rows[(parts[0], parts[1], parts[2], parts[3])] = float(parts[4])
    return rows

with open(out_tsv, "w") as out:
    out.write("sample\tchrom\tstart\tend\tname\tdepth_18h\tdepth_full\tfull_over_18h\n")
    for sample in shared:
        h18 = load(h18_ids[sample])
        ful = load(full_ids[sample])
        for region in h18:
            d18 = h18[region]
            df  = ful.get(region, float("nan"))
            ratio = (df / d18) if d18 and d18 == d18 else float("nan")
            out.write("\t".join([
                sample, region[0], region[1], region[2], region[3],
                f"{d18:.2f}", f"{df:.2f}" if df == df else "NA",
                f"{ratio:.2f}" if ratio == ratio else "NA",
            ]) + "\n")
PYEOF
else
    echo "Skipping paired comparison (need at least one BAM in each group)."
fi
echo

# Summary
SUMMARY="$OUTDIR/v6_cohort_summary.txt"
{
    echo "v6 cohort coverage summary"
    echo "=========================="
    echo
    echo "Panel BED:     $PANEL_BED"
    echo "Panel MD5:     $PANEL_MD5"
    echo "Full-run BAMs: ${#FULL_BAMS[@]}"
    for b in "${FULL_BAMS[@]}"; do echo "    $(basename "$b")"; done
    echo "18h BAMs:      ${#H18_BAMS[@]}"
    for b in "${H18_BAMS[@]}"; do echo "    $(basename "$b")"; done
    echo

    if [[ -s "$COHORT_FULL_TSV" ]]; then
        echo "--- Full-run cohort coverage ($COHORT_FULL_TSV) ---"
        echo "Regions with mean cohort coverage < 5x (SV-calling floor):"
        awk -F'\t' 'NR>1 && $(NF-2)+0 < 5 {printf "  %-15s %s:%s-%s\tmean=%s min=%s max=%s\n",
            $4, $1, $2, $3, $(NF-2), $(NF-1), $NF}' "$COHORT_FULL_TSV"
        echo
        echo "Regions with mean cohort coverage >= 15x (on-spec):"
        awk -F'\t' 'NR>1 && $(NF-2)+0 >= 15 {printf "  %-15s %s:%s-%s\tmean=%s\n",
            $4, $1, $2, $3, $(NF-2)}' "$COHORT_FULL_TSV"
        echo
        echo "Mean coverage across all regions, per full-run sample:"
        awk -F'\t' '
            NR==1 {for(i=5;i<=NF-3;i++) headers[i]=$i; ncols=NF-3; next}
            {for(i=5;i<=ncols;i++) {sum[i]+=$i; n[i]++}}
            END {for(i=5;i<=ncols;i++) printf "  %-30s %.2f\n", headers[i], sum[i]/n[i]}
        ' "$COHORT_FULL_TSV"
        echo
    fi

    if [[ -s "$COHORT_18H_TSV" ]]; then
        echo "--- 18h snapshot coverage ($COHORT_18H_TSV) ---"
        echo "Mean coverage across all regions, per 18h sample:"
        awk -F'\t' '
            NR==1 {for(i=5;i<=NF-3;i++) headers[i]=$i; ncols=NF-3; next}
            {for(i=5;i<=ncols;i++) {sum[i]+=$i; n[i]++}}
            END {for(i=5;i<=ncols;i++) printf "  %-30s %.2f\n", headers[i], sum[i]/n[i]}
        ' "$COHORT_18H_TSV"
        echo
    fi

    if [[ -s "$PAIRED_TSV" ]]; then
        echo "--- 18h vs full paired comparison ($PAIRED_TSV) ---"
        echo "Median (full/18h) ratio per sample (on-target panel regions):"
        echo "  (well-behaved AS run: ratio scales ~linearly with extra time;"
        echo "   ratio plateauing well below run-time-ratio suggests AS losing"
        echo "   efficiency or pores dying.)"
        python3 - "$PAIRED_TSV" <<'PYEOF'
import sys, statistics
from collections import defaultdict
ratios = defaultdict(list)
with open(sys.argv[1]) as fh:
    next(fh)  # header
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8: continue
        sample = parts[0]
        ratio = parts[7]
        if ratio == "NA": continue
        try:
            ratios[sample].append(float(ratio))
        except ValueError:
            continue
for s in sorted(ratios):
    rs = ratios[s]
    if not rs: continue
    print(f"  {s:<30} median ratio = {statistics.median(rs):.2f}  "
          f"(n_regions = {len(rs)}, p25 = {statistics.quantiles(rs, n=4)[0]:.2f}, "
          f"p75 = {statistics.quantiles(rs, n=4)[2]:.2f})")
PYEOF
        echo
    fi

    echo "--- Interpretation notes ---"
    echo "  * v6 regions that are NEW (DIS3, TRAF3, PRDM1, ATM, CYLD, H1-4,"
    echo "    MAX, EGR1, LTB, ATR) were NOT targeted by v5 adaptive sampling."
    echo "    Coverage on these regions reflects off-target shotgun depth"
    echo "    (typically 1-2x). This is the baseline against which v6"
    echo "    sequencing will be compared."
    echo "  * v5-retained regions should match v5 enrichment (15-20x target)."
    echo "  * v6 TP53 window is a strict subset of v5 TP53, so coverage is"
    echo "    unchanged from v5."
} | tee "$SUMMARY"

echo
echo "Done. Outputs in $OUTDIR:"
ls -1 "$OUTDIR"/*.tsv "$OUTDIR"/v6_cohort_summary.txt 2>/dev/null | sed 's|^|  |'
