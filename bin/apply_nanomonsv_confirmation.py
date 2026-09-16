#!/usr/bin/env python3
"""
apply_nanomonsv_confirmation.py

Handoff open item 6 (2026-09-15): wire nanomonsv as a confirmation column.

modules/local/nanomonsv.nf has existed since 14 September, unwired, and
nextflow.config carried none of the params it reads. Its own header sets out
why it is a confirmation layer rather than a fifth caller in the SURVIVOR
merge: with the HPRC control panel it is far more specific than the ensemble
(44 PASS calls genome-wide against ~950 ensemble rows on one sample) and less
sensitive (a FISH-confirmed two-read t(11;14) is captured at parse but never
reported). A junction the ensemble reports and nanomonsv confirms is strongly
supported; one it does not report means nothing either way.

Change
------
  nextflow.config            skip_nanomonsv, nanomonsv_control_panel (HPRC
                             CHM13 prefix), nanomonsv_simple_repeat_bed,
                             nanomonsv_min_read_num 2, nanomonsv_min_vaf 0.02,
                             nanomonsv_tol 100
  conf/envs.config           NANOMONSV runs in awgs_sv
  conf/modules.config        NANOMONSV publishes result.txt and result.vcf to
                             t2t/calls/nanomonsv
  subworkflows/local/t2t_track.nf
                             NANOMONSV(t2t_bam_bai) beside the four callers;
                             its result.txt joined into AUGMENT_SV_SUPPORT
                             with remainder: true, so a sample nanomonsv
                             fails on keeps its annotated table; emits
                             nanomonsv_result
  modules/local/augment_sv_support.nf
                             optional path(nanomonsv_result); passes
                             --nanomonsv when present
  bin/augment_sv_support.py  --nanomonsv <result.txt> [--nanomonsv-tol 100]
                             adds support_nanomonsv, filter_nanomonsv,
                             nanomonsv_confirmed (yes / no / blank when
                             nanomonsv did not run). Not gated on the
                             callers column, not folded into support_reads or
                             filter_worst: it is read beside the ensemble,
                             never as part of it.
  bin/merge_translocations.py
                             a merged junction is confirmed if any member
                             is; support_nanomonsv is the cluster maximum,
                             matching the other support_* columns

The match is orientation-agnostic on both breakends within nanomonsv_tol bp,
default 100: nanomonsv reports the precise junction while SURVIVOR positions
are shifted by up to survivor_max_dist, so the confirmation tolerance is
wider than the per-caller support tolerance (25 bp).

Idempotent; sentinel `nanomonsv_confirmation_v1`; --check; --revert.
"""

import argparse
import os
import shutil
import sys

SENTINEL = "nanomonsv_confirmation_v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

EDITS = []

# ---------------------------------------------------------------------------
EDITS.append(("nextflow.config", [(
    "    support_tol   = 25\n",
    "    support_tol   = 25\n"
    "\n"
    "    // ===== nanomonsv confirmation layer (nanomonsv_confirmation_v1) =====\n"
    "    // Runs beside the four-caller ensemble, not inside the SURVIVOR merge;\n"
    "    // see modules/local/nanomonsv.nf for the measurements behind that.\n"
    "    // The HPRC control panel is what gives it its specificity: V(D)J\n"
    "    // recombination and somatic hypermutation are in every normal B cell\n"
    "    // and the panel subtracts them.\n"
    "    skip_nanomonsv              = false\n"
    "    nanomonsv_control_panel     = '/goast/hemat_data/references/nanomonsv/hprc_year1_data_freeze_nanopore_guppy6_minimap2_2_24_merge_control_CHM13'\n"
    "    nanomonsv_simple_repeat_bed = '/goast/hemat_data/references/nanomonsv/human_chm13v2.0_simpleRepeat.bed.gz'\n"
    "    // Below the tool's defaults (3 reads, 0.05): adaptive sampling puts\n"
    "    // real breakpoints at single-digit support.\n"
    "    nanomonsv_min_read_num      = 2\n"
    "    nanomonsv_min_vaf           = 0.02\n"
    "    // Per-breakend tolerance for matching a nanomonsv junction to an\n"
    "    // ensemble row. Wider than support_tol because SURVIVOR positions\n"
    "    // move by up to survivor_max_dist while nanomonsv's are exact.\n"
    "    nanomonsv_tol               = 100\n",
)]))

# ---------------------------------------------------------------------------
EDITS.append(("conf/envs.config", [(
    "    withName: ICHORKARYO                { conda = params.conda_awgs_sv }\n",
    "    withName: ICHORKARYO                { conda = params.conda_awgs_sv }\n"
    "    withName: NANOMONSV                 { conda = params.conda_awgs_sv }   // nanomonsv_confirmation_v1\n",
)]))

# ---------------------------------------------------------------------------
EDITS.append(("conf/modules.config", [(
    "    withName: SEVERUS {\n",
    "    withName: NANOMONSV {\n"
    "        // nanomonsv_confirmation_v1. The parse step leaves several GB of\n"
    "        // per-read intermediates in the task directory; only the results\n"
    "        // are published.\n"
    "        publishDir = [\n"
    "            path: { \"${params.outdir}/t2t/calls/nanomonsv\" },\n"
    "            mode: 'copy',\n"
    "            pattern: '*.nanomonsv.*'\n"
    "        ]\n"
    "    }\n"
    "    withName: SEVERUS {\n",
)]))

# ---------------------------------------------------------------------------
EDITS.append(("subworkflows/local/t2t_track.nf", [
    (
        "include { SAVANA               } from '../../modules/local/savana'\n",
        "include { SAVANA               } from '../../modules/local/savana'\n"
        "include { NANOMONSV            } from '../../modules/local/nanomonsv'   // nanomonsv_confirmation_v1\n",
    ),
    (
        "        SNIFFLES(t2t_bam_bai)\n"
        "        CUTESV(t2t_bam_bai)\n"
        "        SEVERUS(t2t_bam_bai)\n"
        "        SAVANA(t2t_bam_bai)\n",
        "        SNIFFLES(t2t_bam_bai)\n"
        "        CUTESV(t2t_bam_bai)\n"
        "        SEVERUS(t2t_bam_bai)\n"
        "        SAVANA(t2t_bam_bai)\n"
        "\n"
        "        // Confirmation layer, not a fifth vote: nanomonsv runs beside the\n"
        "        // ensemble and its junctions are matched onto the annotated table\n"
        "        // downstream. It never enters the SURVIVOR merge, so SUPP_VEC and\n"
        "        // the callers column keep their meaning.\n"
        "        nanomonsv_ch = Channel.empty()\n"
        "        if (!params.skip_nanomonsv) {\n"
        "            NANOMONSV(t2t_bam_bai)\n"
        "            nanomonsv_ch = NANOMONSV.out.result\n"
        "        }\n",
    ),
    (
        "            ch_augment_in = ANNOTATE_MM_TRANSLOCATIONS.out.tsv\n"
        "                .join(SNIFFLES.out.vcf, by: 0)\n"
        "                .join(CUTESV.out.vcf,   by: 0)\n"
        "                .join(SEVERUS.out.vcf,  by: 0)\n"
        "                .join(SAVANA.out.vcf,   by: 0)\n"
        "                .map { meta, annotated, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf, sa_vcf ->\n"
        "                    tuple(meta, annotated, sn_vcf, cu_vcf, sv_vcf, sa_vcf)\n"
        "                }\n",
        "            // The nanomonsv result is joined with remainder: true on meta.id,\n"
        "            // so a sample nanomonsv fails on (or a run with it skipped)\n"
        "            // still gets its annotated table, with the confirmation columns\n"
        "            // blank. A plain join would drop the sample silently.\n"
        "            ch_augment_in = ANNOTATE_MM_TRANSLOCATIONS.out.tsv\n"
        "                .join(SNIFFLES.out.vcf, by: 0)\n"
        "                .join(CUTESV.out.vcf,   by: 0)\n"
        "                .join(SEVERUS.out.vcf,  by: 0)\n"
        "                .join(SAVANA.out.vcf,   by: 0)\n"
        "                .map { meta, annotated, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf, sa_vcf ->\n"
        "                    tuple(meta.id, meta, annotated, sn_vcf, cu_vcf, sv_vcf, sa_vcf)\n"
        "                }\n"
        "                .join(nanomonsv_ch.map { meta, res -> tuple(meta.id, res) }, remainder: true)\n"
        "                .filter { it[1] != null }\n"
        "                .map { _id, meta, annotated, sn_vcf, cu_vcf, sv_vcf, sa_vcf, nano ->\n"
        "                    tuple(meta, annotated, sn_vcf, cu_vcf, sv_vcf, sa_vcf, nano ?: [])\n"
        "                }\n",
    ),
    (
        "    merged_vcf       = params.skip_sv_calling     ? Channel.empty() : SURVIVOR_MERGE.out.merged_vcf\n",
        "    merged_vcf       = params.skip_sv_calling     ? Channel.empty() : SURVIVOR_MERGE.out.merged_vcf\n"
        "    nanomonsv_result = (params.skip_sv_calling || params.skip_nanomonsv) ? Channel.empty() : NANOMONSV.out.result\n",
    ),
]))

# ---------------------------------------------------------------------------
EDITS.append(("modules/local/augment_sv_support.nf", [
    (
        "          path(severus_vcf),\n"
        "          path(savana_vcf)\n",
        "          path(severus_vcf),\n"
        "          path(savana_vcf),\n"
        "          path(nanomonsv_result)   // nanomonsv_confirmation_v1; [] when absent\n",
    ),
    (
        "    script:\n"
        "    \"\"\"\n"
        "    augment_sv_support.py \\\\\n"
        "        --annotated annotated_in.tsv \\\\\n"
        "        --sniffles  ${sniffles_vcf} \\\\\n"
        "        --cutesv    ${cutesv_vcf} \\\\\n"
        "        --severus   ${severus_vcf} \\\\\n"
        "        --savana    ${savana_vcf} \\\\\n"
        "        --output    ${meta.id}.mm_annotated.tsv \\\\\n"
        "        --tol       ${params.support_tol}\n",
        "    script:\n"
        "    def nano_arg = nanomonsv_result ? \"--nanomonsv ${nanomonsv_result} --nanomonsv-tol ${params.nanomonsv_tol}\" : ''\n"
        "    \"\"\"\n"
        "    augment_sv_support.py \\\\\n"
        "        --annotated annotated_in.tsv \\\\\n"
        "        --sniffles  ${sniffles_vcf} \\\\\n"
        "        --cutesv    ${cutesv_vcf} \\\\\n"
        "        --severus   ${severus_vcf} \\\\\n"
        "        --savana    ${savana_vcf} \\\\\n"
        "        --output    ${meta.id}.mm_annotated.tsv \\\\\n"
        "        --tol       ${params.support_tol} \\\\\n"
        "        ${nano_arg}\n",
    ),
]))

# ---------------------------------------------------------------------------
EDITS.append(("bin/augment_sv_support.py", [
    (
        "def best_match(ca, pa, cb, pb, records, tol):\n",
        "def load_nanomonsv(path):\n"
        "    \"\"\"nanomonsv_confirmation_v1: (canon ends, supporting reads, Is_Filter)\n"
        "    for every row of a nanomonsv result.txt. Columns are read by header\n"
        "    name; Dir_1/Dir_2 are not used because the ensemble match is\n"
        "    orientation-agnostic. Filtered rows are kept so a junction nanomonsv\n"
        "    saw and rejected can be told apart from one it never reported.\"\"\"\n"
        "    out = []\n"
        "    if not path or not os.path.isfile(path):\n"
        "        return out\n"
        "    with open_any(path) as fh:\n"
        "        reader = csv.DictReader(fh, delimiter=\"\\t\")\n"
        "        for row in reader:\n"
        "            try:\n"
        "                c1, p1 = row[\"Chr_1\"], int(row[\"Pos_1\"])\n"
        "                c2, p2 = row[\"Chr_2\"], int(row[\"Pos_2\"])\n"
        "            except (KeyError, TypeError, ValueError):\n"
        "                continue\n"
        "            sup = row.get(\"Supporting_Read_Num_Tumor\", \"\")\n"
        "            sup = int(sup) if str(sup).isdigit() else None\n"
        "            if sup is None:\n"
        "                continue\n"
        "            filt = (row.get(\"Is_Filter\") or \"\").strip() or \"PASS\"\n"
        "            ca, pa, cb, pb = canon(c1, p1, c2, p2)\n"
        "            out.append((ca, pa, cb, pb, sup, filt))\n"
        "    return out\n"
        "\n"
        "\n"
        "def best_match(ca, pa, cb, pb, records, tol):\n",
    ),
    (
        "    ap.add_argument(\"--output\", required=True)\n"
        "    ap.add_argument(\"--tol\", type=int, default=25, help=\"bp tolerance per breakpoint (default 25).\")\n"
        "    args = ap.parse_args()\n",
        "    ap.add_argument(\"--output\", required=True)\n"
        "    ap.add_argument(\"--tol\", type=int, default=25, help=\"bp tolerance per breakpoint (default 25).\")\n"
        "    # nanomonsv_confirmation_v1\n"
        "    ap.add_argument(\"--nanomonsv\", default=None,\n"
        "                    help=\"nanomonsv result.txt; adds support_nanomonsv, \"\n"
        "                         \"filter_nanomonsv and nanomonsv_confirmed. Read beside \"\n"
        "                         \"the ensemble, never folded into support_reads.\")\n"
        "    ap.add_argument(\"--nanomonsv-tol\", type=int, default=100,\n"
        "                    help=\"bp tolerance per breakend for the confirmation match \"\n"
        "                         \"(default 100; SURVIVOR shifts positions, nanomonsv does not).\")\n"
        "    args = ap.parse_args()\n",
    ),
    (
        "    sav = load_caller(args.savana, \"TUMOUR_READ_SUPPORT\")\n"
        "    sys.stderr.write(\n"
        "        f\"loaded support records: sniffles={len(snf)} cutesv={len(cut)} \"\n"
        "        f\"severus={len(sev)} savana={len(sav)}\\n\")\n",
        "    sav = load_caller(args.savana, \"TUMOUR_READ_SUPPORT\")\n"
        "    nano = load_nanomonsv(args.nanomonsv)\n"
        "    nano_run = bool(args.nanomonsv)\n"
        "    sys.stderr.write(\n"
        "        f\"loaded support records: sniffles={len(snf)} cutesv={len(cut)} \"\n"
        "        f\"severus={len(sev)} savana={len(sav)} nanomonsv={len(nano)}\"\n"
        "        f\"{'' if nano_run else ' (not run)'}\\n\")\n",
    ),
    (
        "                \"filter_savana\",\n"
        "                \"filter_worst\"]\n",
        "                \"filter_savana\",\n"
        "                \"filter_worst\",\n"
        "                \"support_nanomonsv\", \"filter_nanomonsv\", \"nanomonsv_confirmed\"]\n",
    ),
    (
        "        seen = [f for f in (sf, cf, vf) if f]\n"
        "        nonpass = sorted({f for f in seen if f != \"PASS\"})\n"
        "        r[\"filter_worst\"] = \";\".join(nonpass) if nonpass else (\"PASS\" if seen else \"\")\n",
        "        seen = [f for f in (sf, cf, vf) if f]\n"
        "        nonpass = sorted({f for f in seen if f != \"PASS\"})\n"
        "        r[\"filter_worst\"] = \";\".join(nonpass) if nonpass else (\"PASS\" if seen else \"\")\n"
        "        # nanomonsv_confirmation_v1. Not gated on the callers column: the\n"
        "        # point is an independent method seeing the same junction. 'yes'\n"
        "        # needs a PASS nanomonsv junction within tolerance on both ends;\n"
        "        # a filtered nanomonsv hit is 'no' with its filter recorded; blank\n"
        "        # means nanomonsv did not run for this sample.\n"
        "        if nano_run:\n"
        "            ns, nf = best_match(ca, pa, cb, pb, nano, args.nanomonsv_tol)\n"
        "            r[\"support_nanomonsv\"] = str(ns) if ns is not None else \"\"\n"
        "            r[\"filter_nanomonsv\"] = nf or \"\"\n"
        "            r[\"nanomonsv_confirmed\"] = \"yes\" if (ns is not None and nf == \"PASS\") else \"no\"\n"
        "        else:\n"
        "            r[\"support_nanomonsv\"] = r[\"filter_nanomonsv\"] = r[\"nanomonsv_confirmed\"] = \"\"\n",
    ),
]))

# ---------------------------------------------------------------------------
EDITS.append(("bin/merge_translocations.py", [
    (
        "    sup_sev = _max_col(\"support_severus\")\n",
        "    sup_sev = _max_col(\"support_severus\")\n"
        "    # nanomonsv_confirmation_v1: confirmed if any member is; the cluster\n"
        "    # maximum for support, matching the other support_* columns.\n"
        "    sup_nano = _max_col(\"support_nanomonsv\")\n"
        "    nano_values = [(m.get(\"nanomonsv_confirmed\") or \"\").strip() for m in members]\n"
        "    nano_confirmed = (\"yes\" if \"yes\" in nano_values\n"
        "                      else (\"no\" if \"no\" in nano_values else \"\"))\n",
    ),
    (
        "        \"support_severus\": sup_sev,\n"
        "        \"n_merged\": str(len(members)),\n",
        "        \"support_severus\": sup_sev,\n"
        "        \"support_nanomonsv\": sup_nano,\n"
        "        \"nanomonsv_confirmed\": nano_confirmed,\n"
        "        \"n_merged\": str(len(members)),\n",
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
