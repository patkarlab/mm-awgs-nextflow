#!/usr/bin/env bash
#
# run_v7_batches.sh
#
# Run several v7 batches through the pipeline one at a time, verifying each
# before starting the next.
#
# Sequential by design. Each batch stages BAMs into work/ and then copies
# published outputs into results/, so peak usage is roughly twice the final
# size of a run. The 20260713 batch alone finished at 241G, meaning ~500G in
# flight. Two of these concurrently would fill the volume.
#
# Usage:
#   tools/run_v7_batches.sh samplesheets/samplesheet_v7_20260713.csv \
#                           samplesheets/samplesheet_v7_20260721.csv
#
#   MIN_FREE_GB=600 tools/run_v7_batches.sh <sheets...>
#   REAP=1           tools/run_v7_batches.sh <sheets...>   # drop results after verify
#   DRY_RUN=1        tools/run_v7_batches.sh <sheets...>   # plan only
#
# REAP deletes the results tree once the batch verifies, keeping the report
# bundle and its zip. Use it when running the full set: five batches of
# results will not fit alongside each other.

set -euo pipefail

MIN_FREE_GB="${MIN_FREE_GB:-600}"
REAP="${REAP:-0}"
DRY_RUN="${DRY_RUN:-0}"
PROFILE="${PROFILE:-docker}"

if [[ $# -eq 0 ]]; then
    echo "Usage: run_v7_batches.sh <samplesheet.csv> [more.csv ...]" >&2
    exit 1
fi

source /home/hemat/anaconda3/etc/profile.d/conda.sh
conda activate awgs_sv

free_gb() {
    df -BG --output=avail /goast | tail -1 | tr -dc '0-9'
}

# Verify every sheet before launching anything, so a typo in the last one is
# not discovered four hours into the first.
for sheet in "$@"; do
    if [[ ! -f "$sheet" ]]; then
        echo "ERROR: samplesheet not found: $sheet" >&2
        exit 1
    fi
    while IFS=, read -r sample bam rest; do
        [[ "$sample" == "sample_id" || -z "$sample" ]] && continue
        if [[ ! -e "$bam" ]]; then
            echo "ERROR: $sheet: input not found for ${sample}: $bam" >&2
            exit 1
        fi
    done < "$sheet"
    n=$(( $(wc -l < "$sheet") - 1 ))
    echo "OK  $sheet  (${n} samples, all inputs present)"
done

echo ""
echo "Batches   : $#"
echo "Profile   : $PROFILE"
echo "Min free  : ${MIN_FREE_GB}G before each batch"
echo "Reap      : $([[ "$REAP" == "1" ]] && echo "yes, results removed after verify" || echo "no")"
echo "Free now  : $(free_gb)G"
echo ""

if [[ "$DRY_RUN" == "1" ]]; then
    echo "Dry run: nothing launched."
    exit 0
fi

overall=0

for sheet in "$@"; do
    tag=$(basename "$sheet" .csv | sed 's/^samplesheet_//')
    outdir="results_${tag}"
    bundle="report_${outdir}"
    stamp=$(date +%Y%m%d_%H%M%S)
    log="nohup_${tag}_${stamp}"

    echo "============================================================"
    echo "Batch ${tag}   started $(date +%H:%M:%S)"
    echo "============================================================"

    avail=$(free_gb)
    if [[ "$avail" -lt "$MIN_FREE_GB" ]]; then
        echo "STOP: ${avail}G free, need ${MIN_FREE_GB}G." >&2
        echo "      Reap earlier batches or raise MIN_FREE_GB deliberately." >&2
        overall=1
        break
    fi
    echo "  free space : ${avail}G"

    # work/ is cleared between batches. -resume is unreliable on this install
    # because of the docker.userEmulation deprecation, so a stale cache buys
    # nothing and costs hundreds of gigabytes.
    rm -rf work .nextflow

    set +e
    nextflow run main.nf -profile "$PROFILE" \
        --sample_sheet "$sheet" \
        --outdir "$outdir" \
        > "${log}.out" 2> "${log}.err"
    rc=$?
    set -e

    if [[ $rc -ne 0 ]]; then
        echo "  FAILED (exit ${rc}). Tail of ${log}.err:" >&2
        tail -20 "${log}.err" >&2
        echo "  Stopping; later batches not started." >&2
        overall=1
        break
    fi
    echo "  pipeline   : complete"

    # Verification. A green pipeline has repeatedly produced a report that
    # looked finished and was not, so each batch is checked rather than
    # assumed.
    ok=1

    if [[ -d "${outdir}/igv" ]]; then
        python3 bin/check_igv_pages.py "$outdir" > "${log}.igv" 2>&1 || true
        if grep -q "broken\|missing" "${log}.igv"; then
            totals=$(grep '^Totals:' "${log}.igv" || true)
            echo "  igv pages  : ${totals:-see ${log}.igv}"
            if echo "$totals" | grep -qv '0 broken, 0 missing'; then ok=0; fi
        fi
    else
        echo "  igv pages  : none produced" >&2
    fi

    for s in $(awk -F, 'NR>1 && $1!="" {print $1}' "$sheet"); do
        f="${bundle}/${s}/snv/${s}_somaticseq_clinical_final.tsv"
        if [[ -f "$f" ]]; then
            if ! head -1 "$f" | tr '\t' '\n' | grep -qx 'Start'; then
                echo "  ALIAS FAIL : ${s} variant table has no alias columns" >&2
                ok=0
            fi
        fi
    done

    if [[ -f "${bundle}.zip" ]]; then
        echo "  archive    : ${bundle}.zip ($(du -h "${bundle}.zip" | cut -f1))"
    else
        echo "  archive    : MISSING" >&2
        ok=0
    fi

    if [[ "$ok" -ne 1 ]]; then
        echo "  VERIFY FAILED for ${tag}; stopping." >&2
        overall=1
        break
    fi
    echo "  verified   : ok"

    if [[ "$REAP" == "1" ]]; then
        echo "  reaping    : ${outdir} ($(du -sh "$outdir" | cut -f1))"
        rm -rf "$outdir"
    fi

    echo "  finished   $(date +%H:%M:%S), $(free_gb)G free"
    echo ""
done

rm -rf work .nextflow
echo "============================================================"
echo "All batches done. Free: $(free_gb)G"
ls -lh report_*.zip 2>/dev/null || true
exit "$overall"
