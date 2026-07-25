#!/usr/bin/env python3
"""
patch_baf_loh_tab.py

Append the BAF/LOH enhancement include to the end of its tab pane.

The enhancement operates on the rendered table rather than on the markup that
produces it, so this patch does not need to understand the pane's structure:
it inserts one include line before the pane's closing tag.

Safety: the anchor is validated before writing, a timestamped .bak is kept, and
the edit is skipped if already applied. --dry-run writes nothing.

Usage:
    python3 bin/patch_baf_loh_tab.py --dry-run
    python3 bin/patch_baf_loh_tab.py
"""

import argparse
import os
import re
import shutil
import sys
from datetime import datetime

INCLUDE = "  {% include 'baf_loh_enhance.html.j2' %}\n"
SENTINEL = "baf_loh_enhance.html.j2"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", default=".")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    base = os.path.join(args.repo, "bin", "dashboard_builder", "templates")
    tab = os.path.join(base, "baf_loh_tab.html.j2")
    enhance = os.path.join(base, "baf_loh_enhance.html.j2")

    for path in (tab, enhance):
        if not os.path.isfile(path):
            print("ERROR: not found: %s" % path, file=sys.stderr)
            return 2

    with open(tab, encoding="utf-8") as handle:
        text = handle.read()

    if SENTINEL in text:
        print("SKIP  baf_loh_tab: already patched")
        return 0

    # Insert before the pane's final closing div, which is the last one in the
    # file. Matching the last occurrence keeps the include inside the pane.
    index = text.rstrip().rfind("</div>")
    if index == -1:
        print("ERROR: no closing </div> found in baf_loh_tab.html.j2",
              file=sys.stderr)
        print("       Nothing was written.", file=sys.stderr)
        return 2

    line_start = text.rfind("\n", 0, index) + 1
    text = text[:line_start] + INCLUDE + text[line_start:]

    print("EDIT  baf_loh_tab: add the enhancement include")
    if args.dry_run:
        print("\nDry run: nothing written.")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy2(tab, "%s.bak_baf_%s" % (tab, stamp))
    with open(tab, "w", encoding="utf-8") as handle:
        handle.write(text)
    print("\nApplied. Backup written alongside.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
