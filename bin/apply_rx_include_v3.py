#!/usr/bin/env python3
"""
apply_rx_include_v3.py  (rx_include_v3)

The Include boxes on the genome-wide and per-chromosome figure cards did not
reach the Reporting tab. initCNVGallery() wires them, and it looks for the
pane by the id the CNV tab had before dashboard_mm_v3 merged it into
"Copy number & allelic state" (#tab-cnv -> #tab-copy-number), so it has been
returning early since that merge. It now binds to the copy-number pane, with
the document as fallback, so a future rename cannot switch it off again.

One file: bin/dashboard_builder/templates/sample_report.html.j2. Backup
<file>.rx3.bak.

Usage (from the repo root):
  python3 bin/apply_rx_include_v3.py            apply
  python3 bin/apply_rx_include_v3.py --check
  python3 bin/apply_rx_include_v3.py --revert
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "rx_include_v3"
BAK_SUFFIX = ".rx3.bak"
TARGET = "bin/dashboard_builder/templates/sample_report.html.j2"

EDITS = [
    (
        "    const cnvTab = document.getElementById(\"tab-cnv\");\n"
        "    if (!cnvTab || !window.tspipeReporting) return;\n",
        "    // rx_include_v3: the pane was #tab-cnv before dashboard_mm_v3 merged it\n"
        "    // into #tab-copy-number; the old lookup returned null and the figure\n"
        "    // cards' Include boxes were never wired. Fall back to the document.\n"
        "    const cnvTab = document.getElementById(\"tab-copy-number\")\n"
        "                || document.getElementById(\"tab-cnv\") || document;\n"
        "    if (!window.tspipeReporting) return;\n",
    ),
]


def state(root):
    p = root / TARGET
    if not p.is_file():
        return "MISSING FILE"
    text = p.read_text()
    if SENTINEL in text:
        return "applied"
    bad = [a for a, _ in EDITS if text.count(a) != 1]
    if bad:
        return "ANCHOR MISMATCH"
    return "ready"


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1].strip())
    ap.add_argument("--root", default=".")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true")
    g.add_argument("--revert", action="store_true")
    args = ap.parse_args()
    root = Path(args.root).resolve()
    p = root / TARGET
    bak = p.with_name(p.name + BAK_SUFFIX)
    if args.check:
        print("  %-52s %s" % (TARGET, state(root)))
        return
    if args.revert:
        if bak.is_file():
            shutil.copy2(bak, p); bak.unlink()
            print("  %-52s restored from %s" % (TARGET, BAK_SUFFIX))
        else:
            print("  %-52s no %s, skipped" % (TARGET, BAK_SUFFIX))
        return
    s = state(root)
    if s == "applied":
        print("Already applied; nothing to do."); return
    if s != "ready":
        sys.exit("Refusing to apply: %s %s" % (TARGET, s))
    text = p.read_text()
    shutil.copy2(p, bak)
    for a, r in EDITS:
        text = text.replace(a, r, 1)
    assert SENTINEL in text
    p.write_text(text)
    print("  %-52s patched (backup: %s)" % (TARGET, bak.name))


if __name__ == "__main__":
    main()
