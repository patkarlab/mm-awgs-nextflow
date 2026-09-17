#!/usr/bin/env python3
"""
apply_bundle_cn_guard.py  (bundle_cn_kind_v1)

Require the ichorkaryo karyotype JSON and the ICHORCNA_PLOT figures in the
report-bundle publish guard.

Why: the workflow requires the `cnv` kind, which is ichorCNA's own params.txt.
That is one process upstream of what the dashboard reads from
hg38/copy_number/: <id>.ichorkaryo.json (ICHORKARYO: ISCN, ploidy class,
trisomies, arm events) and <id>.cn_genome.png / <id>.cn_grid.png
(ICHORCNA_PLOT). The poll computed `cn` but nothing required it, and the
figures were not checked at all, so either could land after the bundle was
built, exactly as the genebaf files did on the 14-sample run.

Applies on top of bundle_allelic_kind_v1 (its sentinel must be present).

Edits, all string-anchored:

  bin/check_published_outputs.sh   `cn` documented as required-capable; new
                                   kind `cnfig` = <id>.cn_genome.png and
                                   <id>.cn_grid.png under hg38/copy_number/
  workflows/mm_awgs.nf             `cn` required when !skip_ichorcna &&
                                   !skip_ichorkaryo; `cnfig` additionally
                                   when !skip_cn_plot (the gates in
                                   subworkflows/local/hg38_track.nf)
  bin/build_report_bundle.sh       BUNDLE_STRICT audit checks the JSON and
                                   both figures in <bundle>/<sample>/copy_number/

Usage (from the repo root):
  python3 bin/apply_bundle_cn_guard.py            apply
  python3 bin/apply_bundle_cn_guard.py --check    report state, change nothing
  python3 bin/apply_bundle_cn_guard.py --revert   restore the .bak copies

Backups are written as <file>.cn.bak so they do not overwrite the .bak set
left by the allelic patch.
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "bundle_cn_kind_v1"
PREREQ = "bundle_allelic_kind_v1"
BAK_SUFFIX = ".cn.bak"

EDITS = [
    (
        "bin/check_published_outputs.sh",
        [
            (
                "#   required_kinds  comma-separated subset of: snv,sv,cnv,qc,cn,allelic,baf\n",
                "#   required_kinds  comma-separated subset of: snv,sv,cnv,qc,cn,cnfig,allelic,baf\n",
            ),
            (
                "#   cn    <id>*.ichorkaryo.json                              (ICHORKARYO)\n",
                "#   cn    <id>*.ichorkaryo.json                              (ICHORKARYO)\n"
                "#   cnfig <id>.cn_genome.png and <id>.cn_grid.png\n"
                "#         under hg38/copy_number/                          (ICHORCNA_PLOT)\n"
                "#         bundle_cn_kind_v1: cnv only proves ichorCNA's params.txt\n"
                "#         exists; the JSON and figures are two processes later.\n",
            ),
            (
                "GENEBAF_LIST=$(list_files -path '*hg38/allelic*' -name '*.genebaf*')\n",
                "GENEBAF_LIST=$(list_files -path '*hg38/allelic*' -name '*.genebaf*')\n"
                "CNFIG_LIST=$(list_files -path '*hg38/copy_number*' -name '*.cn_*.png')\n",
            ),
            (
                "KINDS=(snv sv cnv qc cn allelic)\n",
                "KINDS=(snv sv cnv qc cn cnfig allelic)\n",
            ),
            (
                "      cn)  has \"$CN_LIST\" \"$id\" && ok=1 || ok=0 ;;\n",
                "      cn)  has \"$CN_LIST\" \"$id\" && ok=1 || ok=0 ;;\n"
                "      cnfig)\n"
                "        ok=1\n"
                "        for f in cn_genome.png cn_grid.png; do\n"
                "          has \"$CNFIG_LIST\" \"${id}.${f}\" || ok=0\n"
                "        done ;;\n",
            ),
        ],
    ),
    (
        "workflows/mm_awgs.nf",
        [
            (
                "            if (!params.skip_genebaf) required_kinds << 'allelic'\n",
                "            if (!params.skip_genebaf) required_kinds << 'allelic'\n"
                "            // bundle_cn_kind_v1: the karyotype JSON and copy-number figures\n"
                "            // are two processes downstream of the params.txt that 'cnv'\n"
                "            // checks; gates mirror subworkflows/local/hg38_track.nf.\n"
                "            if (!params.skip_ichorcna && !params.skip_ichorkaryo) {\n"
                "                required_kinds << 'cn'\n"
                "                if (!params.skip_cn_plot) required_kinds << 'cnfig'\n"
                "            }\n",
            ),
        ],
    ),
    (
        "bin/build_report_bundle.sh",
        [
            (
                "# translocation tables, allelic = the four genebaf files\n"
                "# (bundle_allelic_kind_v1), baf = the cohort screen table.\n",
                "# translocation tables, allelic = the four genebaf files\n"
                "# (bundle_allelic_kind_v1), cn = the ichorkaryo JSON, cnfig = the genome\n"
                "# and grid copy-number figures (bundle_cn_kind_v1), baf = the cohort\n"
                "# screen table.\n",
            ),
            (
                "        [[ -s \"$BUNDLE/$s/allelic/$f\" ]] || gaps+=(\"$s:allelic/$f\")\n"
                "      done\n"
                "    fi\n",
                "        [[ -s \"$BUNDLE/$s/allelic/$f\" ]] || gaps+=(\"$s:allelic/$f\")\n"
                "      done\n"
                "    fi\n"
                "    if [[ \"$require\" == *,cn,* ]]; then\n"
                "      f=\"${s}.ichorkaryo.json\"\n"
                "      [[ -s \"$BUNDLE/$s/copy_number/$f\" ]] || gaps+=(\"$s:copy_number/$f\")\n"
                "    fi\n"
                "    if [[ \"$require\" == *,cnfig,* ]]; then\n"
                "      for f in \"${s}.cn_genome.png\" \"${s}.cn_grid.png\"; do\n"
                "        [[ -s \"$BUNDLE/$s/copy_number/$f\" ]] || gaps+=(\"$s:copy_number/$f\")\n"
                "      done\n"
                "    fi\n",
            ),
        ],
    ),
]


def state(root):
    out = []
    for rel, edits in EDITS:
        p = root / rel
        if not p.is_file():
            out.append((rel, "MISSING FILE"))
            continue
        text = p.read_text()
        if SENTINEL in text:
            out.append((rel, "applied"))
            continue
        if PREREQ not in text:
            out.append((rel, "PREREQUISITE MISSING: apply bundle_allelic_kind_v1 first"))
            continue
        bad = [a for a, _ in edits if text.count(a) != 1]
        if bad:
            out.append((rel, "ANCHOR MISMATCH: %d of %d anchors not unique" % (len(bad), len(edits))))
            for a in bad:
                sys.stderr.write("    count=%d  %r\n" % (text.count(a), a.splitlines()[0]))
            continue
        out.append((rel, "ready"))
    return out


def apply(root):
    st = state(root)
    if any(s not in ("ready", "applied") for _, s in st):
        for rel, s in st:
            print("  %-36s %s" % (rel, s))
        sys.exit("Refusing to apply: fix the problems above first.")
    if all(s == "applied" for _, s in st):
        print("Already applied (sentinel present in every file); nothing to do.")
        return
    for rel, edits in EDITS:
        p = root / rel
        text = p.read_text()
        if SENTINEL in text:
            print("  %-36s already applied, skipped" % rel)
            continue
        bak = p.with_name(p.name + BAK_SUFFIX)
        shutil.copy2(p, bak)
        for anchor, repl in edits:
            text = text.replace(anchor, repl, 1)
        assert SENTINEL in text, "sentinel absent after patching %s" % rel
        p.write_text(text)
        print("  %-36s patched (backup: %s)" % (rel, bak.name))


def revert(root):
    for rel, _ in EDITS:
        p = root / rel
        bak = p.with_name(p.name + BAK_SUFFIX)
        if not bak.is_file():
            print("  %-36s no %s, skipped" % (rel, BAK_SUFFIX))
            continue
        shutil.copy2(bak, p)
        bak.unlink()
        print("  %-36s restored from %s" % (rel, BAK_SUFFIX))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1].strip())
    ap.add_argument("--root", default=".", help="repository root (default: cwd)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true", help="report state only")
    g.add_argument("--revert", action="store_true", help="restore backup copies")
    args = ap.parse_args()
    root = Path(args.root).resolve()

    if args.check:
        for rel, s in state(root):
            print("  %-36s %s" % (rel, s))
        return
    if args.revert:
        revert(root)
        return
    apply(root)


if __name__ == "__main__":
    main()
