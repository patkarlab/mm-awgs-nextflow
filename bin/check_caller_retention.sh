#!/usr/bin/env bash
#
# check_caller_retention.sh
#
# Per sample, count merged junctions in translocations.tsv that carry each
# caller, and the total. Run once before collapse_bnd_mates_v1 and once after
# (on the re-run results), and compare. The expectation after the collapse:
# total junctions equal or fewer (less positional repair needed), 'severus'
# counts equal or higher, never lower. A drop in the severus column means the
# kept mate is in the wrong orientation for SURVIVOR to match against
# Sniffles/CuteSV, and the collapse should be reverted before committing.
#
# Usage:
#   bin/check_caller_retention.sh results_v7 > /tmp/retention_before.tsv
#   ... apply, re-run ...
#   bin/check_caller_retention.sh results_v7 > /tmp/retention_after.tsv
#   paste /tmp/retention_before.tsv /tmp/retention_after.tsv

set -euo pipefail
RESULTS="${1:?Usage: check_caller_retention.sh <results_dir>}"

printf 'sample\tjunctions\tsniffles\tcutesv\tseverus\tsavana\tnanomonsv_yes\n'
for f in $(find "$RESULTS" -name '*.translocations.tsv' | sort); do
  s=$(basename "$f" .translocations.tsv)
  awk -F'\t' -v s="$s" '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    { n++
      c = tolower($h["callers"])
      if (c ~ /sniffles/) sn++
      if (c ~ /cutesv/)   cu++
      if (c ~ /severus/)  se++
      if (c ~ /savana/)   sa++
      if (("nanomonsv_confirmed" in h) && $h["nanomonsv_confirmed"]=="yes") na++
    }
    END { printf "%s\t%d\t%d\t%d\t%d\t%d\t%d\n", s, n, sn, cu, se, sa, na }' "$f"
done
