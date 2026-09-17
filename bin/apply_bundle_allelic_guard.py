#!/usr/bin/env python3
"""
apply_bundle_allelic_guard.py  (bundle_allelic_kind_v1)

Add the `allelic` artefact kind (the four genebaf outputs under
hg38/allelic/) to the report-bundle publish guard.

Why: REPORT_BUNDLE polls bin/check_published_outputs.sh for the *required*
artefact kinds and builds the bundle as soon as those are present. genebaf was
never a kind, so on the 14-sample run the bundle was assembled at 06:07 while
publishDir landed the genebaf files at 06:11. build_report_bundle.sh already
collects them; it was simply given nothing to collect. The dashboard then
rendered the copy-number tab without allelic state or the inline
allele-fraction plots.

Three edits, all string-anchored:

  bin/check_published_outputs.sh   new kind `allelic` = <id>.genebaf.tsv,
                                   .genebaf.json, .genebaf.bins.tsv,
                                   .genebaf.sites.tsv under hg38/allelic/;
                                   column width widened so the matrix aligns
  workflows/mm_awgs.nf             `allelic` required when !params.skip_genebaf
                                   (GENEBAF is gated on that param alone)
  bin/build_report_bundle.sh       BUNDLE_STRICT audit checks the same four
                                   files in <bundle>/<sample>/allelic/

Usage (from the repo root):
  python3 bin/apply_bundle_allelic_guard.py            apply
  python3 bin/apply_bundle_allelic_guard.py --check    report state, change nothing
  python3 bin/apply_bundle_allelic_guard.py --revert   restore the .bak copies

Each file is written once, after a .bak copy is taken. The sentinel string
`bundle_allelic_kind_v1` refuses a second application.
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "bundle_allelic_kind_v1"

# (path, [(anchor, replacement), ...])
# Every anchor must occur exactly once in the file before patching.
EDITS = [
    (
        "bin/check_published_outputs.sh",
        [
            (
                "#   required_kinds  comma-separated subset of: snv,sv,cnv,qc,cn,baf\n",
                "#   required_kinds  comma-separated subset of: snv,sv,cnv,qc,cn,allelic,baf\n",
            ),
            (
                "#   cn    <id>*.ichorkaryo.json                              (ICHORKARYO)\n",
                "#   cn    <id>*.ichorkaryo.json                              (ICHORKARYO)\n"
                "#   allelic  <id>.genebaf.{tsv,json,bins.tsv,sites.tsv}\n"
                "#         under hg38/allelic/                              (GENEBAF)\n"
                "#         bundle_allelic_kind_v1: added after the 14-sample bundle\n"
                "#         was built four minutes before genebaf's publishDir copy.\n",
            ),
            (
                "CN_LIST=$(list_files -name '*.ichorkaryo.json')\n",
                "CN_LIST=$(list_files -name '*.ichorkaryo.json')\n"
                "GENEBAF_LIST=$(list_files -path '*hg38/allelic*' -name '*.genebaf*')\n",
            ),
            (
                "KINDS=(snv sv cnv qc cn)\n",
                "KINDS=(snv sv cnv qc cn allelic)\n",
            ),
            (
                "      cn)  has \"$CN_LIST\" \"$id\" && ok=1 || ok=0 ;;\n",
                "      cn)  has \"$CN_LIST\" \"$id\" && ok=1 || ok=0 ;;\n"
                "      allelic)\n"
                "        # All four genebaf files, matched on the exact basename so\n"
                "        # <id>.genebaf.tsv does not satisfy the check for .bins.tsv.\n"
                "        ok=1\n"
                "        for suf in .tsv .json .bins.tsv .sites.tsv; do\n"
                "          has \"$GENEBAF_LIST\" \"${id}.genebaf${suf}\" || ok=0\n"
                "        done ;;\n",
            ),
            # Column width: "allelic*" is eight characters and broke the matrix.
            (
                "  if is_required \"$k\"; then printf ' %5s' \"$k*\"; else printf ' %5s' \"$k\"; fi\n",
                "  if is_required \"$k\"; then printf ' %8s' \"$k*\"; else printf ' %8s' \"$k\"; fi\n",
            ),
            (
                "      printf ' %5s' \"+\"\n",
                "      printf ' %8s' \"+\"\n",
            ),
            (
                "      printf ' %5s' \"-\"\n",
                "      printf ' %8s' \"-\"\n",
            ),
        ],
    ),
    (
        "workflows/mm_awgs.nf",
        [
            (
                "            if (!params.skip_ichorcna) required_kinds << 'cnv'\n",
                "            if (!params.skip_ichorcna) required_kinds << 'cnv'\n"
                "            // bundle_allelic_kind_v1: genebaf outputs. GENEBAF reads the\n"
                "            // BAM directly and is gated on skip_genebaf alone.\n"
                "            if (!params.skip_genebaf) required_kinds << 'allelic'\n",
            ),
        ],
    ),
    (
        "bin/build_report_bundle.sh",
        [
            (
                "# to audit (default snv,sv): snv = both aliased variant tables, sv = both\n"
                "# translocation tables, baf = the cohort screen table.\n",
                "# to audit (default snv,sv): snv = both aliased variant tables, sv = both\n"
                "# translocation tables, allelic = the four genebaf files\n"
                "# (bundle_allelic_kind_v1), baf = the cohort screen table.\n",
            ),
            (
                "      for f in \"${s}.mm_annotated.tsv\" \"${s}.translocations.tsv\"; do\n"
                "        [[ -f \"$BUNDLE/$s/translocations/$f\" ]] || gaps+=(\"$s:translocations/$f\")\n"
                "      done\n"
                "    fi\n",
                "      for f in \"${s}.mm_annotated.tsv\" \"${s}.translocations.tsv\"; do\n"
                "        [[ -f \"$BUNDLE/$s/translocations/$f\" ]] || gaps+=(\"$s:translocations/$f\")\n"
                "      done\n"
                "    fi\n"
                "    if [[ \"$require\" == *,allelic,* ]]; then\n"
                "      for f in \"${s}.genebaf.tsv\" \"${s}.genebaf.json\" \\\n"
                "               \"${s}.genebaf.bins.tsv\" \"${s}.genebaf.sites.tsv\"; do\n"
                "        [[ -s \"$BUNDLE/$s/allelic/$f\" ]] || gaps+=(\"$s:allelic/$f\")\n"
                "      done\n"
                "    fi\n",
            ),
        ],
    ),
]


def state(root):
    """Return (path, status) for each target: 'applied', 'ready', or a problem."""
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
        bak = p.with_suffix(p.suffix + ".bak")
        shutil.copy2(p, bak)
        for anchor, repl in edits:
            text = text.replace(anchor, repl, 1)
        assert SENTINEL in text, "sentinel absent after patching %s" % rel
        p.write_text(text)
        print("  %-36s patched (backup: %s)" % (rel, bak.name))


def revert(root):
    for rel, _ in EDITS:
        p = root / rel
        bak = p.with_suffix(p.suffix + ".bak")
        if not bak.is_file():
            print("  %-36s no .bak, skipped" % rel)
            continue
        shutil.copy2(bak, p)
        bak.unlink()
        print("  %-36s restored from .bak" % rel)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1].strip())
    ap.add_argument("--root", default=".", help="repository root (default: cwd)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true", help="report state only")
    g.add_argument("--revert", action="store_true", help="restore .bak copies")
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
