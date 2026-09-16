#!/usr/bin/env python3
"""
apply_genebaf_scope_span.py

Handoff open items 3 and 4 (2026-09-15), in bin/genebaf.py, modules/local/
genebaf.nf and nextflow.config.

Item 4 — BAF scope (decided 15 Sep, not implemented until now)
--------------------------------------------------------------
Depth and log2 stay on every panel window. Allelic state is measured only at
the windows given by --baf-regions (a BED; a panel window is in scope when it
overlaps any region in it), and never at IGH/IGK/IGL: V(D)J recombination and
somatic hypermutation make allele fractions there clonal, not germline. The
exclusion is by window name (substrings, default IGH,IGK,IGL) and applies
whether or not --baf-regions is given, so a panel that forgets to list the
Ig loci still cannot pass them through.

Out-of-scope windows are still in every output with their depth, and carry
allelic_state "not_assessed" plus a baf_scope column that says why
("outside_baf_regions" or "excluded_ig"). "insufficient" keeps its meaning:
assessed, too few hets. Sites are loaded only for in-scope windows, so the
VCF scan and the pileup both shrink with the scope.

Item 3 — covered span on arm verdicts
-------------------------------------
The arm LOH test pools panel windows on an arm. On 1p that was NRAS and
TENT5C, 3 Mb apart and both inside a 35 Mb deletion, so a verdict labelled
"1p" described one segment of it. Each arm verdict now records what the
pooled windows actually cover: span_start, span_end, covered_bp, arm_bp,
covered_fraction (windows / arm), span_fraction (extent / arm). With BAF
scoped to eight windows there will be more such arms, not fewer, and the
reader needs to see it. --min-arm-windows still gates the test; the span
does not change any verdict.

Asset
-----
assets/aWGS_PCN_v7_baf_regions_hg38.bed, delivered alongside: the eight wide
windows (CCND1, CCND2, CCND3, MAFA, MAFB, TP53+TNFSF12, MYC, WWOX/MAF) copied
from assets/aWGS_PCN_v7_hg38.bed. Overlap-based, so it also scopes a v4/v5/v6
sample whose panel names differ.

Idempotent; sentinel `genebaf_scope_span_v1`; --check; --revert.
"""

import argparse
import os
import shutil
import sys

SENTINEL = "genebaf_scope_span_v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ASSET = "assets/aWGS_PCN_v7_baf_regions_hg38.bed"

EDITS = []

# ---------------------------------------------------------------------------
# bin/genebaf.py
# ---------------------------------------------------------------------------
EDITS.append(("bin/genebaf.py", [
    # 1. Scope helpers, after load_bed().
    (
        "def load_sites_in_windows(vcf_path, windows, af_field, min_af, max_af):\n",
        "def load_regions(path):\n"
        "    \"\"\"{chrom: [(start, end), ...]} from a BED; name column ignored.\"\"\"\n"
        "    regions = defaultdict(list)\n"
        "    with open(path, encoding=\"utf-8\") as handle:\n"
        "        for line in handle:\n"
        "            if not line.strip() or line.startswith((\"#\", \"track\", \"browser\")):\n"
        "                continue\n"
        "            parts = line.rstrip(\"\\n\").split(\"\\t\")\n"
        "            if len(parts) < 3:\n"
        "                continue\n"
        "            chrom = canonical(parts[0])\n"
        "            if chrom is None:\n"
        "                continue\n"
        "            try:\n"
        "                regions[chrom].append((int(parts[1]), int(parts[2])))\n"
        "            except ValueError:\n"
        "                continue\n"
        "    return regions\n"
        "\n"
        "\n"
        "def baf_scope_of(window, baf_regions, exclude_names):\n"
        "    \"\"\"genebaf_scope_span_v1: why a window is or is not assessed for BAF.\n"
        "\n"
        "    Returns 'assessed', 'excluded_ig' or 'outside_baf_regions'. The Ig\n"
        "    exclusion is applied first and unconditionally: V(D)J recombination\n"
        "    and somatic hypermutation make allele fractions at IGH/IGK/IGL clonal\n"
        "    rather than germline, so they carry no allelic-state information\n"
        "    whatever the region list says. With no region list every other window\n"
        "    is assessed, which is the pre-scope behaviour.\n"
        "    \"\"\"\n"
        "    chrom, start, end, name = window\n"
        "    upper = name.upper()\n"
        "    if any(token and token in upper for token in exclude_names):\n"
        "        return \"excluded_ig\"\n"
        "    if baf_regions is None:\n"
        "        return \"assessed\"\n"
        "    for r_start, r_end in baf_regions.get(chrom, ()):\n"
        "        if start < r_end and r_start < end:\n"
        "            return \"assessed\"\n"
        "    return \"outside_baf_regions\"\n"
        "\n"
        "\n"
        "def load_sites_in_windows(vcf_path, windows, af_field, min_af, max_af):\n",
    ),
    # 2. arm_loh: pool only assessed windows; record covered span.
    (
        "    pooled = {}\n"
        "    for row in rows:\n"
        "        arm = arm_of(row[\"chrom\"], row[\"start\"], row[\"end\"], arm_bounds)\n"
        "        if arm is None or row[\"n_usable\"] == 0:\n"
        "            continue\n"
        "        entry = pooled.setdefault(arm, {\"chrom\": row[\"chrom\"], \"n_windows\": 0,\n"
        "                                        \"n_het\": 0, \"n_usable\": 0,\n"
        "                                        \"expected\": [], \"genes\": []})\n"
        "        entry[\"n_windows\"] += 1\n"
        "        entry[\"n_het\"] += row[\"n_het\"]\n"
        "        entry[\"n_usable\"] += row[\"n_usable\"]\n"
        "        entry[\"expected\"].append(row[\"expected_het_fraction\"] * row[\"n_usable\"])\n"
        "        entry[\"genes\"].append(row[\"gene\"])\n",
        "    pooled = {}\n"
        "    for row in rows:\n"
        "        # genebaf_scope_span_v1: a window outside the BAF scope has no\n"
        "        # sites and must not count toward n_windows either.\n"
        "        if row.get(\"baf_scope\", \"assessed\") != \"assessed\":\n"
        "            continue\n"
        "        arm = arm_of(row[\"chrom\"], row[\"start\"], row[\"end\"], arm_bounds)\n"
        "        if arm is None or row[\"n_usable\"] == 0:\n"
        "            continue\n"
        "        entry = pooled.setdefault(arm, {\"chrom\": row[\"chrom\"], \"n_windows\": 0,\n"
        "                                        \"n_het\": 0, \"n_usable\": 0,\n"
        "                                        \"expected\": [], \"genes\": [],\n"
        "                                        \"span_start\": row[\"start\"],\n"
        "                                        \"span_end\": row[\"end\"],\n"
        "                                        \"covered_bp\": 0})\n"
        "        entry[\"n_windows\"] += 1\n"
        "        entry[\"n_het\"] += row[\"n_het\"]\n"
        "        entry[\"n_usable\"] += row[\"n_usable\"]\n"
        "        entry[\"expected\"].append(row[\"expected_het_fraction\"] * row[\"n_usable\"])\n"
        "        entry[\"genes\"].append(row[\"gene\"])\n"
        "        entry[\"span_start\"] = min(entry[\"span_start\"], row[\"start\"])\n"
        "        entry[\"span_end\"] = max(entry[\"span_end\"], row[\"end\"])\n"
        "        entry[\"covered_bp\"] += max(0, row[\"end\"] - row[\"start\"])\n",
    ),
    (
        "        out[arm] = {\n"
        "            \"chrom\": entry[\"chrom\"],\n"
        "            \"n_windows\": entry[\"n_windows\"],\n"
        "            \"genes\": entry[\"genes\"],\n",
        "        # genebaf_scope_span_v1: what the pooled windows actually cover, so a\n"
        "        # verdict labelled '1p' can be read as 'two windows 3 Mb apart on 1p'.\n"
        "        # Arm length is p: chrom start to centromere, q: centromere to end.\n"
        "        chrom_start, cen_start, cen_end, chrom_end = arm_bounds[entry[\"chrom\"]]\n"
        "        arm_bp = (cen_start - chrom_start) if arm.endswith(\"p\") else (chrom_end - cen_end)\n"
        "        span_bp = entry[\"span_end\"] - entry[\"span_start\"]\n"
        "        out[arm] = {\n"
        "            \"chrom\": entry[\"chrom\"],\n"
        "            \"n_windows\": entry[\"n_windows\"],\n"
        "            \"genes\": entry[\"genes\"],\n"
        "            \"span_start\": entry[\"span_start\"],\n"
        "            \"span_end\": entry[\"span_end\"],\n"
        "            \"covered_bp\": entry[\"covered_bp\"],\n"
        "            \"arm_bp\": arm_bp,\n"
        "            \"covered_fraction\": (round(entry[\"covered_bp\"] / arm_bp, 4) if arm_bp > 0 else None),\n"
        "            \"span_fraction\": (round(span_bp / arm_bp, 4) if arm_bp > 0 else None),\n",
    ),
    # 3. CLI.
    (
        "    parser.add_argument(\"--max-other-fraction\", type=float, default=0.2,\n"
        "                        help=\"drop sites where third-allele reads exceed this \"\n"
        "                             \"fraction of depth [0.2]\")\n"
        "    args = parser.parse_args()\n",
        "    parser.add_argument(\"--max-other-fraction\", type=float, default=0.2,\n"
        "                        help=\"drop sites where third-allele reads exceed this \"\n"
        "                             \"fraction of depth [0.2]\")\n"
        "    # genebaf_scope_span_v1\n"
        "    parser.add_argument(\"--baf-regions\", default=None,\n"
        "                        help=\"BED of regions at which allelic state is measured. \"\n"
        "                             \"A panel window is assessed when it overlaps any \"\n"
        "                             \"region; the rest keep depth and log2 only. Without \"\n"
        "                             \"it every non-Ig window is assessed.\")\n"
        "    parser.add_argument(\"--baf-exclude\", default=\"IGH,IGK,IGL\",\n"
        "                        help=\"comma-separated window-name substrings never \"\n"
        "                             \"assessed for allelic state, whatever --baf-regions \"\n"
        "                             \"says [IGH,IGK,IGL]\")\n"
        "    args = parser.parse_args()\n",
    ),
    # 4. Scope the windows before loading sites.
    (
        "    sites = load_sites_in_windows(args.sites_vcf, windows, args.af_field,\n"
        "                                  args.min_af, args.max_af)\n"
        "    total_sites = sum(len(v) for v in sites.values())\n"
        "    sys.stderr.write(f\"{total_sites} common SNPs inside panel windows\\n\")\n",
        "    # genebaf_scope_span_v1: decide the BAF scope first, and load sites only\n"
        "    # for windows that will be assessed.\n"
        "    baf_regions = load_regions(args.baf_regions) if args.baf_regions else None\n"
        "    exclude_names = [t.strip().upper() for t in args.baf_exclude.split(\",\") if t.strip()]\n"
        "    scope = {w: baf_scope_of(w, baf_regions, exclude_names) for w in windows}\n"
        "    baf_windows = [w for w in windows if scope[w] == \"assessed\"]\n"
        "    sys.stderr.write(f\"{len(baf_windows)} windows in BAF scope\"\n"
        "                     f\" ({sum(1 for s in scope.values() if s == 'excluded_ig')} Ig excluded,\"\n"
        "                     f\" {sum(1 for s in scope.values() if s == 'outside_baf_regions')}\"\n"
        "                     f\" outside --baf-regions)\\n\")\n"
        "\n"
        "    sites = load_sites_in_windows(args.sites_vcf, baf_windows, args.af_field,\n"
        "                                  args.min_af, args.max_af)\n"
        "    total_sites = sum(len(v) for v in sites.values())\n"
        "    sys.stderr.write(f\"{total_sites} common SNPs inside BAF-scope windows\\n\")\n",
    ),
    # 5. Skip the pileup for out-of-scope windows.
    (
        "        site_list = sites.get(window, [])\n"
        "        counts = tally_window(bam, chrom, start, end, site_list,\n"
        "                              args.min_mapq, args.min_baseq)\n",
        "        site_list = sites.get(window, []) if scope[window] == \"assessed\" else []\n"
        "        counts = (tally_window(bam, chrom, start, end, site_list,\n"
        "                               args.min_mapq, args.min_baseq)\n"
        "                  if site_list else {})\n",
    ),
    (
        "        state, stat = classify(n_het, n_hom, expected_het_fraction, devs, in_band,\n"
        "                               args.min_hets, args.loh_ratio)\n"
        "\n"
        "        rows.append({\n"
        "            \"gene\": name, \"chrom\": chrom, \"start\": start, \"end\": end,\n",
        "        if scope[window] == \"assessed\":\n"
        "            state, stat = classify(n_het, n_hom, expected_het_fraction, devs, in_band,\n"
        "                                   args.min_hets, args.loh_ratio)\n"
        "        else:\n"
        "            state, stat = \"not_assessed\", None\n"
        "\n"
        "        rows.append({\n"
        "            \"gene\": name, \"chrom\": chrom, \"start\": start, \"end\": end,\n"
        "            \"baf_scope\": scope[window],\n",
    ),
    # 6. Output columns.
    (
        "            handle.write(\"arm\\tchrom\\tn_windows\\tn_usable\\tn_het\\t\"\n"
        "                         \"expected_het_fraction\\tobserved_het_fraction\\tratio\\t\"\n"
        "                         \"p_value\\tverdict\\tgenes\\n\")\n"
        "            for arm, v in arms.items():\n"
        "                handle.write(\"\\t\".join(str(x) for x in (\n"
        "                    arm, v[\"chrom\"], v[\"n_windows\"], v[\"n_usable\"], v[\"n_het\"],\n"
        "                    v[\"expected_het_fraction\"], v[\"observed_het_fraction\"],\n"
        "                    \"\" if v[\"ratio\"] is None else v[\"ratio\"],\n"
        "                    \"\" if v[\"p_value\"] is None else v[\"p_value\"],\n"
        "                    v[\"verdict\"], \",\".join(v[\"genes\"]))) + \"\\n\")\n",
        "            handle.write(\"arm\\tchrom\\tn_windows\\tn_usable\\tn_het\\t\"\n"
        "                         \"expected_het_fraction\\tobserved_het_fraction\\tratio\\t\"\n"
        "                         \"p_value\\tverdict\\tspan_start\\tspan_end\\tcovered_bp\\t\"\n"
        "                         \"arm_bp\\tcovered_fraction\\tspan_fraction\\tgenes\\n\")\n"
        "            for arm, v in arms.items():\n"
        "                handle.write(\"\\t\".join(str(x) for x in (\n"
        "                    arm, v[\"chrom\"], v[\"n_windows\"], v[\"n_usable\"], v[\"n_het\"],\n"
        "                    v[\"expected_het_fraction\"], v[\"observed_het_fraction\"],\n"
        "                    \"\" if v[\"ratio\"] is None else v[\"ratio\"],\n"
        "                    \"\" if v[\"p_value\"] is None else v[\"p_value\"],\n"
        "                    v[\"verdict\"],\n"
        "                    v[\"span_start\"], v[\"span_end\"], v[\"covered_bp\"], v[\"arm_bp\"],\n"
        "                    \"\" if v[\"covered_fraction\"] is None else v[\"covered_fraction\"],\n"
        "                    \"\" if v[\"span_fraction\"] is None else v[\"span_fraction\"],\n"
        "                    \",\".join(v[\"genes\"]))) + \"\\n\")\n",
    ),
    (
        "        cols = [\"gene\", \"chrom\", \"start\", \"end\", \"mean_depth\", \"n_sites\", \"n_usable\",\n",
        "        cols = [\"gene\", \"chrom\", \"start\", \"end\", \"baf_scope\", \"mean_depth\", \"n_sites\", \"n_usable\",\n",
    ),
    (
        "        \"states\": {s: sum(1 for r in rows if r[\"allelic_state\"] == s)\n"
        "                   for s in (\"balanced\", \"imbalanced\", \"insufficient\")},\n",
        "        \"states\": {s: sum(1 for r in rows if r[\"allelic_state\"] == s)\n"
        "                   for s in (\"balanced\", \"imbalanced\", \"insufficient\", \"not_assessed\")},\n"
        "        \"baf_scope\": {\n"
        "            \"regions\": args.baf_regions,\n"
        "            \"exclude\": exclude_names,\n"
        "            \"assessed_windows\": [w[3] for w in baf_windows],\n"
        "            \"n_assessed\": len(baf_windows),\n"
        "            \"n_excluded_ig\": sum(1 for s in scope.values() if s == \"excluded_ig\"),\n"
        "            \"n_outside_regions\": sum(1 for s in scope.values() if s == \"outside_baf_regions\"),\n"
        "        },\n",
    ),
]))

# ---------------------------------------------------------------------------
# modules/local/genebaf.nf
# ---------------------------------------------------------------------------
EDITS.append(("modules/local/genebaf.nf", [
    (
        "    \"\"\"\n"
        "    set -euo pipefail\n"
        "\n"
        "    genebaf.py \\\\\n",
        "    // genebaf_scope_span_v1: allelic state at the eight wide windows only;\n"
        "    // depth on all. IGH/IGK/IGL are excluded inside genebaf.py regardless.\n"
        "    def baf_regions_arg = params.genebaf_baf_regions ? \"--baf-regions ${params.genebaf_baf_regions}\" : ''\n"
        "    \"\"\"\n"
        "    set -euo pipefail\n"
        "\n"
        "    genebaf.py \\\\\n",
    ),
    (
        "        --max-arm-p ${params.genebaf_max_arm_p}\n",
        "        --max-arm-p ${params.genebaf_max_arm_p} \\\\\n"
        "        ${baf_regions_arg}\n",
    ),
]))

# ---------------------------------------------------------------------------
# nextflow.config
# ---------------------------------------------------------------------------
EDITS.append(("nextflow.config", [
    (
        "    genebaf_max_arm_p       = 0.001\n",
        "    genebaf_max_arm_p       = 0.001\n"
        "    // genebaf_scope_span_v1: BED of windows at which allelic state is\n"
        "    // measured (the eight wide translocation/17p windows). Depth and log2\n"
        "    // are still computed on every panel window. null assesses every\n"
        "    // non-Ig window.\n"
        "    genebaf_baf_regions     = \"${projectDir}/assets/aWGS_PCN_v7_baf_regions_hg38.bed\"\n",
    ),
]))


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def check():
    problems = 0
    if os.path.isfile(os.path.join(REPO, ASSET)):
        print(f"PRESENT  {ASSET}")
    else:
        print(f"MISSING  {ASSET} (copy it from the inbox first)")
        problems += 1
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            print(f"MISSING  {rel}"); problems += 1; continue
        text = read(path)
        if SENTINEL in text:
            print(f"APPLIED  {rel}"); continue
        bad = [old for old, _ in pairs if text.count(old) != 1]
        if bad:
            print(f"NO-MATCH {rel}: {len(bad)} anchor(s) not found exactly once")
            for old in bad:
                print(f"         anchor starts: {old.strip().splitlines()[0]!r}")
            problems += 1
        else:
            print(f"READY    {rel}")
    return 1 if problems else 0


def apply():
    if not os.path.isfile(os.path.join(REPO, ASSET)):
        sys.exit(f"ERROR: {ASSET} not found; copy it into place first. Nothing changed.")
    staged = []
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {rel} not found; nothing changed.")
        text = read(path)
        if SENTINEL in text:
            print(f"skip     {rel} (already applied)"); continue
        for old, new in pairs:
            n = text.count(old)
            if n != 1:
                sys.exit(f"ERROR: anchor found {n} times in {rel} (expected 1); nothing changed.\n"
                         f"       anchor starts: {old.strip().splitlines()[0]!r}")
            text = text.replace(old, new, 1)
        if SENTINEL not in text:
            sys.exit(f"ERROR: internal: sentinel absent after edit of {rel}")
        staged.append((rel, path, text))
    for rel, path, text in staged:
        shutil.copy2(path, path + ".bak")
        write(path, text)
        print(f"patched  {rel} (backup: {rel}.bak)")


def revert():
    for rel, _ in EDITS:
        path = os.path.join(REPO, rel)
        if os.path.isfile(path + ".bak"):
            shutil.move(path + ".bak", path); print(f"restored {rel}")
        else:
            print(f"no .bak  {rel}")


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true")
    g.add_argument("--revert", action="store_true")
    args = ap.parse_args()
    if args.check:
        sys.exit(check())
    if args.revert:
        revert(); return
    apply()


if __name__ == "__main__":
    main()
