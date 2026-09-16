#!/usr/bin/env python3
"""
apply_bundle_guard.py

Idempotent patch for the report-bundle publish race (open items 9 and 10 of
the 2026-09-15 handoff).

Defect. REPORT_BUNDLE assembles the bundle by scanning the published results
tree, and is ordered after the producing processes through a channel of
completion signals. That ordering is against task completion, not against
publishDir: Nextflow copies published files asynchronously after a task
finishes, so the bundle can start while copies are still landing. A run that
was resumed makes this worse, because every cached task re-publishes at once
when the run starts. The symptom is a subset of samples missing tables the
pipeline demonstrably produced (7 of 15 without variant tables while
FILTER_V6_REPORT reads 15/15), and the cohort baf_loh/ directory missing
although the screen completed. BAF_CN_PLOTS was also never in the ordering
signal at all, so its figures could be missing even without the race.

Fix, in four files:

  nextflow.config          params.bundle_publish_wait_s, params.bundle_strict
  subworkflows/local/hg38_track.nf
                           emit baf_cn_figures so the bundle can wait on it
  workflows/mm_awgs.nf     pass the sample ID list and required kinds into
                           REPORT_BUNDLE; add BAF_CN_PLOTS to the signal
  modules/local/report_bundle.nf
                           poll bin/check_published_outputs.sh until every
                           required artefact is present or the deadline
                           passes; export BUNDLE_STRICT
  bin/build_report_bundle.sh
                           with BUNDLE_STRICT=1, fail after assembly if any
                           sample lacks its SNV tables or SV tables

bin/check_published_outputs.sh is a new file delivered alongside this script;
it must be in bin/ before the patched pipeline runs.

Conventions: string-match anchors, .bak written before every edit, sentinel
string refuses a double apply, explicit file list, --check and --revert.

Usage:
  python3 bin/apply_bundle_guard.py            # apply
  python3 bin/apply_bundle_guard.py --check    # report state, change nothing
  python3 bin/apply_bundle_guard.py --revert   # restore from .bak
"""

import argparse
import os
import shutil
import sys

SENTINEL = "bundle_guard_v1"

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# ---------------------------------------------------------------------------
# Edits. Each is (path, [(old, new), ...]). Every `old` must occur exactly
# once in the file; every `new` carries the sentinel somewhere in the file so
# a second apply is refused.
# ---------------------------------------------------------------------------

EDITS = []

# 1. nextflow.config -------------------------------------------------------
EDITS.append((
    "nextflow.config",
    [(
        "    report_bundle_name    = null    // null means derive from outdir\n",
        "    report_bundle_name    = null    // null means derive from outdir\n"
        "    // bundle_guard_v1: publishDir copies land asynchronously after task\n"
        "    // completion, so REPORT_BUNDLE waits until every required artefact is\n"
        "    // visible in the published tree before scanning it.\n"
        "    bundle_publish_wait_s = 900     // seconds to wait for published outputs\n"
        "    bundle_strict         = true    // fail the bundle if any sample lacks required tables\n",
    )],
))

# 2. hg38_track.nf -----------------------------------------------------------
EDITS.append((
    "subworkflows/local/hg38_track.nf",
    [
        (
            "    baf_screen_ch = Channel.empty()\n"
            "    if (!params.skip_clair3_phased && !params.skip_baf_loh) {\n",
            "    baf_screen_ch  = Channel.empty()\n"
            "    baf_figures_ch = Channel.empty()   // bundle_guard_v1\n"
            "    if (!params.skip_clair3_phased && !params.skip_baf_loh) {\n",
        ),
        (
            "            BAF_CN_PLOTS(\n"
            "                BAF_LOH_SCREEN.out.screen,\n"
            "                BAF_LOH_SCREEN.out.sample_map,\n"
            "                clair3_ch,\n"
            "                ichor_ch,\n"
            "                file(params.cohort_bed_hg38, checkIfExists: true)\n"
            "            )\n",
            "            BAF_CN_PLOTS(\n"
            "                BAF_LOH_SCREEN.out.screen,\n"
            "                BAF_LOH_SCREEN.out.sample_map,\n"
            "                clair3_ch,\n"
            "                ichor_ch,\n"
            "                file(params.cohort_bed_hg38, checkIfExists: true)\n"
            "            )\n"
            "            baf_figures_ch = BAF_CN_PLOTS.out.figures\n",
        ),
        (
            "    baf_loh_screen           = baf_screen_ch\n",
            "    baf_loh_screen           = baf_screen_ch\n"
            "    baf_cn_figures           = baf_figures_ch\n",
        ),
    ],
))

# 3. mm_awgs.nf ------------------------------------------------------------
EDITS.append((
    "workflows/mm_awgs.nf",
    [
        (
            "                .mix(HG38_TRACK.out.baf_loh_screen.map  { it -> 'ok' })\n"
            "        }\n",
            "                .mix(HG38_TRACK.out.baf_loh_screen.map  { it -> 'ok' })\n"
            "                // bundle_guard_v1: the per-sample BAF figures were never in\n"
            "                // the signal, so the bundle could start before they existed.\n"
            "                .mix(HG38_TRACK.out.baf_cn_figures.map  { it -> 'ok' })\n"
            "        }\n",
        ),
        (
            "        REPORT_BUNDLE(\n"
            "            ready.collect().ifEmpty(['none']),\n",
            "        // bundle_guard_v1: the bundle is handed the full sample list and the\n"
            "        // artefact kinds this configuration must have produced, and waits\n"
            "        // until all of them are visible in the published tree. Completion\n"
            "        // signals alone are not enough: publishDir copies land after the\n"
            "        // task completes, not before its outputs are emitted.\n"
            "        sample_ids = PREPARE_INPUT.out.minknow_bams\n"
            "            .map { meta, _input -> meta.id.toString() }\n"
            "            .collect()\n"
            "        required_kinds = []\n"
            "        if (!params.skip_t2t_track) {\n"
            "            required_kinds << 'sv'\n"
            "            if (!params.skip_qc) required_kinds << 'qc'\n"
            "        }\n"
            "        if (!params.skip_hg38_track) {\n"
            "            if (!params.skip_clair3_phased && !params.skip_vep_annotate\n"
            "                    && !params.skip_v6_filter) required_kinds << 'snv'\n"
            "            if (!params.skip_ichorcna) required_kinds << 'cnv'\n"
            "            if (!params.skip_clair3_phased && !params.skip_baf_loh) required_kinds << 'baf'\n"
            "        }\n"
            "\n"
            "        REPORT_BUNDLE(\n"
            "            ready.collect().ifEmpty(['none']),\n"
            "            sample_ids,\n"
            "            required_kinds.join(','),\n",
        ),
    ],
))

# 4. report_bundle.nf ---------------------------------------------------------
EDITS.append((
    "modules/local/report_bundle.nf",
    [
        (
            "    val  ready_signals\n"
            "    path bundle_script\n",
            "    val  ready_signals\n"
            "    // bundle_guard_v1: every sample in the run, and the artefact kinds the\n"
            "    // current configuration must have published for each of them.\n"
            "    val  sample_ids\n"
            "    val  required_kinds\n"
            "    path bundle_script\n",
        ),
        (
            "    bash ${bundle_script} \"${results_dir}\" \"${bundle_name}\"\n",
            "    # bundle_guard_v1. Completion of the producing tasks is already\n"
            "    # guaranteed by ready_signals. What is not guaranteed is that publishDir\n"
            "    # has finished copying their outputs into the results tree, and on a\n"
            "    # resumed run every cached task re-publishes at once. Poll until the\n"
            "    # tree holds every required artefact for every sample, then fail\n"
            "    # loudly rather than ship a bundle with silent gaps.\n"
            "    printf '%s\\\\n' ${sample_ids.join(' ')} > sample_ids.txt\n"
            "    deadline=\\$(( SECONDS + ${params.bundle_publish_wait_s} ))\n"
            "    until check_published_outputs.sh \"${results_dir}\" sample_ids.txt \"${required_kinds}\" \\\\\n"
            "            > publish_check.log 2>&1; do\n"
            "        if (( SECONDS >= deadline )); then\n"
            "            cat publish_check.log >&2\n"
            "            echo \"ERROR: required outputs still absent from ${results_dir} after\" >&2\n"
            "            echo \"       ${params.bundle_publish_wait_s} s. Either publishDir has not\" >&2\n"
            "            echo \"       written them or a producing process did not emit them.\" >&2\n"
            "            echo \"       Re-run bin/check_published_outputs.sh by hand to see which.\" >&2\n"
            "            exit 1\n"
            "        fi\n"
            "        sleep 30\n"
            "    done\n"
            "    cat publish_check.log\n"
            "\n"
            "    export BUNDLE_STRICT=\"${params.bundle_strict ? '1' : '0'}\"\n"
            "    export BUNDLE_REQUIRE=\"${required_kinds}\"\n"
            "\n"
            "    bash ${bundle_script} \"${results_dir}\" \"${bundle_name}\"\n",
        ),
    ],
))

# 5. build_report_bundle.sh ------------------------------------------------
EDITS.append((
    "bin/build_report_bundle.sh",
    [(
        "echo \"\"\n"
        "echo \"Bundle tree: $BUNDLE/\"\n",
        "# ---------------------------------------------------------------------------\n"
        "# Completeness (bundle_guard_v1).\n"
        "#\n"
        "# Every collector above reports a missing source as a warning and moves on,\n"
        "# which is right for optional material but let a bundle with 7 of 15\n"
        "# samples lacking variant tables build, zip and ship as if complete. With\n"
        "# BUNDLE_STRICT=1 the assembled tree is audited against what this run must\n"
        "# contain and the script fails on any gap. BUNDLE_REQUIRE lists the kinds\n"
        "# to audit (default snv,sv): snv = both aliased variant tables, sv = both\n"
        "# translocation tables, baf = the cohort screen table.\n"
        "if [[ \"${BUNDLE_STRICT:-0}\" == \"1\" ]]; then\n"
        "  require=\",${BUNDLE_REQUIRE:-snv,sv},\"\n"
        "  gaps=()\n"
        "  for s in \"${SAMPLES[@]}\"; do\n"
        "    if [[ \"$require\" == *,snv,* ]]; then\n"
        "      for f in \"${s}_somaticseq_clinical_final.tsv\" \"${s}_somaticseq_filtered.tsv\"; do\n"
        "        [[ -s \"$BUNDLE/$s/snv/$f\" ]] || gaps+=(\"$s:snv/$f\")\n"
        "      done\n"
        "    fi\n"
        "    if [[ \"$require\" == *,sv,* ]]; then\n"
        "      for f in \"${s}.mm_annotated.tsv\" \"${s}.translocations.tsv\"; do\n"
        "        [[ -f \"$BUNDLE/$s/translocations/$f\" ]] || gaps+=(\"$s:translocations/$f\")\n"
        "      done\n"
        "    fi\n"
        "  done\n"
        "  if [[ \"$require\" == *,baf,* ]]; then\n"
        "    [[ -s \"$BUNDLE/baf_loh/cohort.baf_screen.tsv\" ]] || gaps+=(\"cohort:baf_loh/cohort.baf_screen.tsv\")\n"
        "  fi\n"
        "  if [[ ${#gaps[@]} -gt 0 ]]; then\n"
        "    echo \"\" >&2\n"
        "    echo \"ERROR: bundle is incomplete (BUNDLE_STRICT=1). Missing:\" >&2\n"
        "    printf '  %s\\n' \"${gaps[@]}\" >&2\n"
        "    echo \"Set BUNDLE_STRICT=0 to build anyway.\" >&2\n"
        "    exit 1\n"
        "  fi\n"
        "  echo \"Completeness check passed for ${#SAMPLES[@]} sample(s) (${BUNDLE_REQUIRE:-snv,sv}).\"\n"
        "fi\n"
        "\n"
        "echo \"\"\n"
        "echo \"Bundle tree: $BUNDLE/\"\n",
    )],
))


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def check():
    """Report per-file state without changing anything. Returns exit code."""
    problems = 0
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            print(f"MISSING  {rel}")
            problems += 1
            continue
        text = read(path)
        if SENTINEL in text:
            print(f"APPLIED  {rel}")
            continue
        bad = [old for old, _ in pairs if text.count(old) != 1]
        if bad:
            print(f"NO-MATCH {rel}: {len(bad)} anchor(s) not found exactly once")
            for old in bad:
                first = old.strip().splitlines()[0]
                print(f"         anchor starts: {first!r}")
            problems += 1
        else:
            print(f"READY    {rel}")
    return 1 if problems else 0


def apply():
    # Validate everything before touching anything, so a bad anchor in file 4
    # does not leave files 1-3 half-patched.
    staged = []
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        if not os.path.isfile(path):
            sys.exit(f"ERROR: {rel} not found; nothing changed.")
        text = read(path)
        if SENTINEL in text:
            print(f"skip     {rel} (already applied)")
            continue
        for old, new in pairs:
            n = text.count(old)
            if n != 1:
                first = old.strip().splitlines()[0]
                sys.exit(f"ERROR: anchor found {n} times in {rel} (expected 1); "
                         f"nothing changed.\n       anchor starts: {first!r}")
            text = text.replace(old, new, 1)
        if SENTINEL not in text:
            sys.exit(f"ERROR: internal: sentinel absent after edit of {rel}")
        staged.append((rel, path, text))

    for rel, path, text in staged:
        shutil.copy2(path, path + ".bak")
        write(path, text)
        print(f"patched  {rel} (backup: {rel}.bak)")

    print("\nNext:")
    print("  cp ~/inbox/from_claude/check_published_outputs.sh bin/ && chmod +x bin/check_published_outputs.sh")
    print("  nextflow run . -profile gandalf -stub-run --sample_sheet <sheet> --outdir stub_out")
    print("  git add nextflow.config subworkflows/local/hg38_track.nf workflows/mm_awgs.nf \\")
    print("          modules/local/report_bundle.nf bin/build_report_bundle.sh \\")
    print("          bin/check_published_outputs.sh bin/apply_bundle_guard.py")


def revert():
    for rel, _ in EDITS:
        path = os.path.join(REPO, rel)
        bak = path + ".bak"
        if os.path.isfile(bak):
            shutil.move(bak, path)
            print(f"restored {rel}")
        else:
            print(f"no .bak  {rel}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true", help="report state; change nothing")
    g.add_argument("--revert", action="store_true", help="restore every file from its .bak")
    args = ap.parse_args()
    if args.check:
        sys.exit(check())
    if args.revert:
        revert()
        return
    apply()


if __name__ == "__main__":
    main()
