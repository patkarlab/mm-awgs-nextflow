#!/usr/bin/env python3
"""
apply_retire_mmkaryo.py

Handoff open item 5 (2026-09-15): retire mmkaryo, mmbaf and the COPY_NUMBER
subworkflow; draw the copy-number figures from ichorCNA on the hg38 track.

Before this change the copy-number tab showed ichorkaryo's tables (ichorCNA,
hg38, purity from flow, ploidy pinned at 2) above figures drawn by mmkaryo
(T2T, its own purity/ploidy fit). Two callers, two references, one tab, and
the figure could contradict the numbers it sat above.

After this change:

  modules/local/ichorcna_plot.nf   NEW process ICHORCNA_PLOT: runs
                                   ichorcna_to_mmplot.py then mmplot.py on
                                   ichorCNA's .cna.seg/.seg.txt, sex and
                                   subtitle from the ichorkaryo JSON, no BAF
                                   panel. Publishes cn_genome / cn_grid /
                                   cn_chr*.png to hg38/copy_number beside
                                   the ichorkaryo JSON.
  subworkflows/local/hg38_track.nf ICHORCNA_PLOT wired after ICHORKARYO;
                                   emits cn_figures.
  workflows/mm_awgs.nf             COPY_NUMBER include, call and bundle
                                   signals removed; cn_figures added to the
                                   bundle signal.
  conf/modules.config              MMKARYO/MMBAF/MMPLOT publishDir blocks
                                   removed; ICHORCNA_PLOT added.
  conf/envs.config                 ICHORCNA_PLOT gets the awgs_sv env.
  nextflow.config                  skip_copy_number retired; skip_cn_plot,
                                   cn_plot_max_cn, cn_plot_chromosomes,
                                   gene_model_hg38 added. The mmkaryo_* and
                                   mmbaf_* params are left in place, unused,
                                   so a saved -params-file still loads.
  bin/build_report_bundle.sh       the four mmkaryo collectors
                                   (karyotype.json, cytobands.tsv,
                                   segments.baf.tsv, arms.tsv) removed; the
                                   PNG copy_all scoped to hg38/copy_number.
  subworkflows/local/copy_number.nf, modules/local/mm{karyo,baf,plot}.nf
                                   RETIRED header prepended; files kept for
                                   one release so history stays legible.

Not touched: bin/mmkaryo.py, bin/mmbaf.py, bin/mmplot.py (mmplot is still
the renderer), bin/dashboard_builder/parsers/copy_number.py (already reads
the ichorkaryo JSON and the cn_*.png names, which do not change).

Prerequisite: modules/local/ichorcna_plot.nf must be in place before the
patched pipeline runs. This script checks for it.

Conventions: string-match anchors, every file validated before any file is
written, .bak per file, sentinel refuses a double apply, --check, --revert.
"""

import argparse
import os
import shutil
import sys

SENTINEL = "retire_mmkaryo_v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NEW_MODULE = "modules/local/ichorcna_plot.nf"

EDITS = []

# ---------------------------------------------------------------------------
# 1. hg38_track.nf
# ---------------------------------------------------------------------------
EDITS.append(("subworkflows/local/hg38_track.nf", [
    (
        "include { ICHORKARYO           } from '../../modules/local/ichorkaryo.nf'\n",
        "include { ICHORKARYO           } from '../../modules/local/ichorkaryo.nf'\n"
        "include { ICHORCNA_PLOT        } from '../../modules/local/ichorcna_plot.nf'   // retire_mmkaryo_v1\n",
    ),
    (
        "                    .map { _id, meta, seg, p, c, w, bed ->\n"
        "                           tuple(meta, seg, p, c, w, bed) }\n"
        "            )\n"
        "        }\n"
        "    }\n",
        "                    .map { _id, meta, seg, p, c, w, bed ->\n"
        "                           tuple(meta, seg, p, c, w, bed) }\n"
        "            )\n"
        "\n"
        "            // Copy-number figures from the same ichorCNA output the\n"
        "            // karyotype JSON was derived from, so the figure and the\n"
        "            // tables on the copy-number tab can never disagree. Replaces\n"
        "            // MMPLOT on mmkaryo (T2T). Joined on meta.id; a sample that\n"
        "            // lacks either ichorCNA file is dropped here, which is the\n"
        "            // right outcome for a figure.\n"
        "            if (!params.skip_cn_plot) {\n"
        "                ICHORCNA_PLOT(\n"
        "                    ICHORCNA.out.cna_seg.map { meta, c -> tuple(meta.id, meta, c) }\n"
        "                        .join(ICHORCNA.out.seg.map { meta, s -> tuple(meta.id, s) })\n"
        "                        .join(ICHORKARYO.out.karyotype.map { meta, k -> tuple(meta.id, k) })\n"
        "                        .join(ch_panel_bed.map { meta, bed -> tuple(meta.id, bed) })\n"
        "                        .map { _id, meta, c, s, k, bed -> tuple(meta, c, s, k, bed) },\n"
        "                    file(params.cytoband_txt_hg38, checkIfExists: true),\n"
        "                    params.gene_model_hg38 ? file(params.gene_model_hg38, checkIfExists: true)\n"
        "                                           : file('NO_FILE')\n"
        "                )\n"
        "            }\n"
        "        }\n"
        "    }\n",
    ),
    (
        "    ichorkaryo_json          = (params.skip_ichorcna || params.skip_ichorkaryo)\n"
        "                               ? Channel.empty() : ICHORKARYO.out.karyotype\n",
        "    ichorkaryo_json          = (params.skip_ichorcna || params.skip_ichorkaryo)\n"
        "                               ? Channel.empty() : ICHORKARYO.out.karyotype\n"
        "    cn_figures               = (params.skip_ichorcna || params.skip_ichorkaryo || params.skip_cn_plot)\n"
        "                               ? Channel.empty() : ICHORCNA_PLOT.out.genome\n",
    ),
]))

# ---------------------------------------------------------------------------
# 2. mm_awgs.nf
# ---------------------------------------------------------------------------
EDITS.append(("workflows/mm_awgs.nf", [
    (
        "include { COPY_NUMBER   } from '../subworkflows/local/copy_number.nf'\n",
        "// COPY_NUMBER (mmkaryo/mmbaf/mmplot on T2T) retired 2026-09-16; copy-number\n"
        "// figures now come from ICHORCNA_PLOT inside HG38_TRACK. retire_mmkaryo_v1\n",
    ),
    (
        "    // 4a. Allele-specific copy number from the off-target reads.\n"
        "    //\n"
        "    // Runs on the T2T track: the panel BED, the centromere table and the\n"
        "    // cytoband file are all CHM13v2.0 coordinates. Distinct from ichorCNA on\n"
        "    // the hg38 track, which is a tumour-fraction estimator that happens to\n"
        "    // emit copy number; this produces segmented allele-specific states, a\n"
        "    // per-segment clonal cell fraction comparable to FISH percentages, and an\n"
        "    // explicit detection limit.\n"
        "    //\n"
        "    // Purity is read from meta.purity (sample sheet). Depth alone cannot\n"
        "    // separate purity from ploidy, and leaving it free put two cohort samples\n"
        "    // at 0.26 and 0.30 against measured flow values of 0.88 and 0.98.\n"
        "    if (!params.skip_copy_number && !params.skip_t2t_track) {\n"
        "        COPY_NUMBER(\n"
        "            T2T_TRACK.out.t2t_bam_bai,\n"
        "            PREPARE_INPUT.out.panel_beds.map { meta, bed_t2t, bed_hg38 ->\n"
        "                tuple(meta, bed_t2t)\n"
        "            },\n"
        "            file(params.t2t_fai,   checkIfExists: true),\n"
        "            file(params.t2t_fasta, checkIfExists: true),\n"
        "            params.t2t_cytobands  ? file(params.t2t_cytobands)  : file('NO_FILE'),\n"
        "            params.t2t_gene_model ? file(params.t2t_gene_model) : file('NO_FILE'),\n"
        "            file(params.g1000_t2t_vcf, checkIfExists: true)\n"
        "        )\n"
        "    }\n"
        "\n",
        "    // 4a. Copy number. Depth-based copy number is ichorCNA on the hg38 track,\n"
        "    // read into a karyotype by ICHORKARYO and drawn by ICHORCNA_PLOT, both\n"
        "    // inside HG38_TRACK. Allelic state at the eight wide panel windows is\n"
        "    // GENEBAF, also inside HG38_TRACK. The former T2T COPY_NUMBER subworkflow\n"
        "    // (mmkaryo/mmbaf/mmplot) is retired: it masked the panel before\n"
        "    // segmenting, which removed focal del(17p) from the data, and its free\n"
        "    // purity/ploidy fit contradicted the flow-purity, diploid-pinned tables.\n"
        "\n",
    ),
    (
        "        // The bundle scans the published tree, so it must not start before the\n"
        "        // copy number figures and tables have been published.\n"
        "        if (!params.skip_copy_number && !params.skip_t2t_track) {\n"
        "            ready = ready\n"
        "                .mix(COPY_NUMBER.out.karyotype.map   { it -> 'ok' })\n"
        "                .mix(COPY_NUMBER.out.plot_genome.map { it -> 'ok' })\n"
        "                .mix(COPY_NUMBER.out.baf.map         { it -> 'ok' })\n"
        "        }\n",
        "",
    ),
    (
        "                .mix(HG38_TRACK.out.baf_cn_figures.map  { it -> 'ok' })\n"
        "        }\n",
        "                .mix(HG38_TRACK.out.baf_cn_figures.map  { it -> 'ok' })\n"
        "                // Copy-number figures and the karyotype JSON are both\n"
        "                // collected by the bundle.\n"
        "                .mix(HG38_TRACK.out.ichorkaryo_json.map { it -> 'ok' })\n"
        "                .mix(HG38_TRACK.out.cn_figures.map      { it -> 'ok' })\n"
        "        }\n",
    ),
]))

# ---------------------------------------------------------------------------
# 3. conf/modules.config
# ---------------------------------------------------------------------------
EDITS.append(("conf/modules.config", [
    (
        "    // ===== Allele-specific copy number (T2T track) =====\n"
        "    // Published under t2t/ because the panel BED, centromere table and\n"
        "    // cytoband file are all CHM13v2.0 coordinates. Kept separate from\n"
        "    // hg38/cnv, which holds ichorCNA: that is a tumour-fraction estimator\n"
        "    // that emits copy number, this is segmented allele-specific copy number.\n"
        "    withName: MMKARYO {\n"
        "        publishDir = [\n"
        "            path: { \"${params.outdir}/t2t/copy_number\" },\n"
        "            mode: 'copy'\n"
        "        ]\n"
        "    }\n"
        "    withName: MMBAF {\n"
        "        publishDir = [\n"
        "            path: { \"${params.outdir}/t2t/copy_number\" },\n"
        "            mode: 'copy'\n"
        "        ]\n"
        "    }\n"
        "    withName: MMPLOT {\n"
        "        publishDir = [\n"
        "            path: { \"${params.outdir}/t2t/copy_number\" },\n"
        "            mode: 'copy',\n"
        "            pattern: '*.png'\n"
        "        ]\n"
        "    }\n"
        "\n",
        "    // MMKARYO / MMBAF / MMPLOT (t2t/copy_number) retired 2026-09-16;\n"
        "    // see ICHORCNA_PLOT under the hg38 track. retire_mmkaryo_v1\n"
        "\n",
    ),
    (
        "    withName: ICHORKARYO {\n"
        "        publishDir = [\n"
        "            path: { \"${params.outdir}/hg38/copy_number\" },\n"
        "            mode: 'copy'\n"
        "        ]\n"
        "    }\n",
        "    withName: ICHORKARYO {\n"
        "        publishDir = [\n"
        "            path: { \"${params.outdir}/hg38/copy_number\" },\n"
        "            mode: 'copy'\n"
        "        ]\n"
        "    }\n"
        "    withName: ICHORCNA_PLOT {\n"
        "        // Beside the ichorkaryo JSON: one directory holds the numbers and\n"
        "        // the figures drawn from them.\n"
        "        publishDir = [\n"
        "            path: { \"${params.outdir}/hg38/copy_number\" },\n"
        "            mode: 'copy',\n"
        "            pattern: '*.{png,tsv}'\n"
        "        ]\n"
        "    }\n",
    ),
]))

# ---------------------------------------------------------------------------
# 4. conf/envs.config
# ---------------------------------------------------------------------------
EDITS.append(("conf/envs.config", [
    (
        "    withName: MMKARYO                   { conda = params.conda_awgs_sv }\n"
        "    withName: MMBAF                     { conda = params.conda_awgs_sv }\n"
        "    withName: MMPLOT                    { conda = params.conda_awgs_sv }\n",
        "    // MMKARYO / MMBAF / MMPLOT retired 2026-09-16 (retire_mmkaryo_v1); a\n"
        "    // selector with no process only produces a warning, but three of\n"
        "    // them on every run is noise.\n"
        "    withName: ICHORCNA_PLOT             { conda = params.conda_awgs_sv }\n",
    ),
]))

# ---------------------------------------------------------------------------
# 5. nextflow.config
# ---------------------------------------------------------------------------
EDITS.append(("nextflow.config", [
    (
        "    skip_copy_number      = false\n",
        "    // skip_copy_number retired 2026-09-16 (retire_mmkaryo_v1): the T2T\n"
        "    // mmkaryo/mmbaf/mmplot subworkflow no longer runs. Copy-number figures\n"
        "    // are ICHORCNA_PLOT on the hg38 track, controlled below. The mmkaryo_*\n"
        "    // and mmbaf_* params that follow are unused and kept only so a saved\n"
        "    // -params-file still loads.\n"
        "    skip_copy_number      = true\n"
        "    skip_cn_plot          = false\n"
        "    cn_plot_max_cn        = 6\n"
        "    cn_plot_chromosomes   = '1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y'\n"
        "    // Gene model in hg38 coordinates for highlight bands on the per-\n"
        "    // chromosome pages. Optional; build with bin/build_mm_gene_model.py\n"
        "    // from the hg38 RefSeq GFF. null draws no gene track.\n"
        "    gene_model_hg38       = null\n",
    ),
]))

# ---------------------------------------------------------------------------
# 6. bin/build_report_bundle.sh
# ---------------------------------------------------------------------------
EDITS.append(("bin/build_report_bundle.sh", [
    (
        "  copy_all   \"$d/copy_number\" \"$s\" -path '*copy_number*' -name '*.png'\n"
        "  copy_first \"$d/copy_number\" \"${s}.ichorkaryo.json\"  \"$s\" -path '*copy_number*' -name '*.ichorkaryo.json'\n"
        "  copy_first \"$d/copy_number\" \"${s}.karyotype.json\"   \"$s\" -path '*copy_number*' -name '*.karyotype.json'\n"
        "  copy_first \"$d/copy_number\" \"${s}.cytobands.tsv\"    \"$s\" -path '*copy_number*' -name '*.cytobands.tsv'\n"
        "  copy_first \"$d/copy_number\" \"${s}.segments.baf.tsv\" \"$s\" -path '*copy_number*' -name '*.segments.baf.tsv'\n"
        "  copy_first \"$d/copy_number\" \"${s}.arms.tsv\"         \"$s\" -path '*copy_number*' -name '*.arms.tsv'\n",
        "  # retire_mmkaryo_v1: figures and JSON both come from hg38/copy_number now.\n"
        "  # The path is scoped so a stale t2t/copy_number directory from a run that\n"
        "  # predates the change cannot leak mmkaryo figures into this tab, and the\n"
        "  # four mmkaryo tables (karyotype.json, cytobands.tsv, segments.baf.tsv,\n"
        "  # arms.tsv) are no longer collected.\n"
        "  copy_all   \"$d/copy_number\" \"$s\" -path '*hg38/copy_number*' -name '*.png'\n"
        "  copy_first \"$d/copy_number\" \"${s}.ichorkaryo.json\"  \"$s\" -path '*hg38/copy_number*' -name '*.ichorkaryo.json'\n",
    ),
]))

# ---------------------------------------------------------------------------
# 7. Retired files: header only.
# ---------------------------------------------------------------------------
RETIRED_HEADER = (
    "// RETIRED 2026-09-16 (retire_mmkaryo_v1). No longer included by any\n"
    "// workflow. Copy-number figures come from ICHORCNA_PLOT on the hg38 track;\n"
    "// see modules/local/ichorcna_plot.nf. Kept for one release so the history\n"
    "// of the copy-number tab stays readable, then delete.\n"
)
for rel, first_line in (
    ("subworkflows/local/copy_number.nf", "//\n// COPY_NUMBER: allele-specific copy number from the off-target reads of the\n"),
    ("modules/local/mmkaryo.nf", "process MMKARYO {\n"),
    ("modules/local/mmbaf.nf", "process MMBAF {\n"),
    ("modules/local/mmplot.nf", "process MMPLOT {\n"),
):
    EDITS.append((rel, [(first_line, RETIRED_HEADER + first_line)]))


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def write(path, text):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)


def check():
    problems = 0
    if not os.path.isfile(os.path.join(REPO, NEW_MODULE)):
        print(f"MISSING  {NEW_MODULE} (copy it from the inbox first)")
        problems += 1
    else:
        print(f"PRESENT  {NEW_MODULE}")
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
                print(f"         anchor starts: {old.strip().splitlines()[0]!r}")
            problems += 1
        else:
            print(f"READY    {rel}")
    return 1 if problems else 0


def apply():
    if not os.path.isfile(os.path.join(REPO, NEW_MODULE)):
        sys.exit(f"ERROR: {NEW_MODULE} not found; copy it into place first. Nothing changed.")
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
    print("\nNext:")
    print("  nextflow run . -profile stub -stub-run --outdir stub_out   # DAG check")
    print("  then resume the v7 session; ICHORCNA_PLOT runs for 15 samples, the rest is cached")


def revert():
    for rel, _ in EDITS:
        path = os.path.join(REPO, rel)
        if os.path.isfile(path + ".bak"):
            shutil.move(path + ".bak", path)
            print(f"restored {rel}")
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
        revert()
        return
    apply()


if __name__ == "__main__":
    main()
