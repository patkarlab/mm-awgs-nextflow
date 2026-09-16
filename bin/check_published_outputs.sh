#!/usr/bin/env bash
#
# check_published_outputs.sh
#
# Verify that the published results tree holds every per-sample artefact the
# report bundle collects, for every sample in the run. Prints a sample x kind
# matrix and exits non-zero if anything is missing.
#
# Two uses:
#
#   1. Diagnostic, run by hand. Answers "is the file absent from the tree, or
#      is the bundle failing to find it?" in one command. For the 7-of-15
#      missing variant tables, a "-" under snv here means the file was never
#      published; a "+" means the bundle ran before publishDir had copied it.
#
#   2. Guard inside REPORT_BUNDLE. publishDir copies run asynchronously after
#      a task completes, and a process that reads the published tree can start
#      before those copies land. The process polls this script until it passes
#      or a deadline expires.
#
# The lookup patterns are the same name-based `find ... | grep -F <id>`
# patterns the bundle uses, so a pass here means the bundle will find the
# same files. Paths under the results directory are not assumed.
#
# Usage:
#   check_published_outputs.sh <results_dir> <ids_file> [required_kinds]
#
#   ids_file        one sample ID per line
#   required_kinds  comma-separated subset of: snv,sv,cnv,qc,cn,baf
#                   (default: snv,sv). Anything not listed is reported but
#                   does not affect the exit status.
#
# Kinds:
#   snv   <id>*.v6_clinical.tsv and <id>*.v6_filtered.tsv   (FILTER_V6_REPORT)
#   sv    <id>.mm_annotated.tsv and <id>.translocations.tsv (T2T track)
#   cnv   ichorCNA params.txt for <id>                       (ICHORCNA)
#   qc    <id>.region_coverage.tsv                           (QC_ONTARGET)
#   cn    <id>*.ichorkaryo.json                              (ICHORKARYO)
#   baf   cohort.baf_screen.tsv, once per run                (BAF_LOH_SCREEN)

set -euo pipefail

RESULTS="${1:?Usage: check_published_outputs.sh <results_dir> <ids_file> [required_kinds]}"
IDS_FILE="${2:?Usage: check_published_outputs.sh <results_dir> <ids_file> [required_kinds]}"
REQUIRED="${3:-snv,sv}"

if [[ ! -d "$RESULTS" ]]; then
  echo "ERROR: results dir not found: $RESULTS" >&2
  exit 2
fi
if [[ ! -f "$IDS_FILE" ]]; then
  echo "ERROR: ids file not found: $IDS_FILE" >&2
  exit 2
fi

mapfile -t IDS < <(grep -v '^[[:space:]]*$' "$IDS_FILE" | sed 's/[[:space:]]*$//' | sort -u)
if [[ ${#IDS[@]} -eq 0 ]]; then
  echo "ERROR: no sample IDs in $IDS_FILE" >&2
  exit 2
fi

# One find per pattern over the whole tree, then per-sample grep against the
# cached list. On a 15-sample results tree this is a few hundred paths.
list_files() {
  # list_files <find-args...>
  find "$RESULTS" -type f "$@" 2>/dev/null || true
}

CLIN_LIST=$(list_files -name '*clinical.tsv')
FILT_LIST=$(list_files -name '*filtered.tsv')
MMANN_LIST=$(list_files -name '*.mm_annotated.tsv')
TRA_LIST=$(list_files -name '*.translocations.tsv')
ICHOR_LIST=$(list_files -path '*ichor*' -name '*params.txt')
QC_LIST=$(list_files -name '*.region_coverage.tsv')
CN_LIST=$(list_files -name '*.ichorkaryo.json')
BAF_LIST=$(list_files -name 'cohort.baf_screen.tsv')

has() {
  # has <list> <id>  -> 0 if any path in <list> contains <id>
  local list="$1" id="$2"
  [[ -n "$list" ]] && grep -qF -- "$id" <<< "$list"
}

is_required() {
  [[ ",${REQUIRED}," == *",$1,"* ]]
}

KINDS=(snv sv cnv qc cn)
missing=()

printf '%-16s' "sample"
for k in "${KINDS[@]}"; do
  if is_required "$k"; then printf ' %5s' "$k*"; else printf ' %5s' "$k"; fi
done
printf '\n'

for id in "${IDS[@]}"; do
  printf '%-16s' "$id"
  for k in "${KINDS[@]}"; do
    case "$k" in
      snv) has "$CLIN_LIST" "$id" && has "$FILT_LIST" "$id" && ok=1 || ok=0 ;;
      sv)  has "$MMANN_LIST" "$id" && has "$TRA_LIST" "$id" && ok=1 || ok=0 ;;
      cnv) has "$ICHOR_LIST" "$id" && ok=1 || ok=0 ;;
      qc)  has "$QC_LIST" "$id" && ok=1 || ok=0 ;;
      cn)  has "$CN_LIST" "$id" && ok=1 || ok=0 ;;
    esac
    if [[ $ok -eq 1 ]]; then
      printf ' %5s' "+"
    else
      printf ' %5s' "-"
      if is_required "$k"; then missing+=("${id}:${k}"); fi
    fi
  done
  printf '\n'
done

# Cohort-level artefact, checked once.
if [[ -n "$BAF_LIST" ]]; then
  echo "baf (cohort)     +   $(head -1 <<< "$BAF_LIST")"
else
  echo "baf (cohort)     -"
  if is_required baf; then missing+=("cohort:baf"); fi
fi

echo ""
echo "Required kinds: ${REQUIRED} (marked *). Samples: ${#IDS[@]}."

if [[ ${#missing[@]} -gt 0 ]]; then
  echo "MISSING (${#missing[@]}): ${missing[*]}" >&2
  exit 1
fi
echo "All required artefacts present."
exit 0
