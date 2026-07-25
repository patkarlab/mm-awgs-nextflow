#!/usr/bin/env python3
"""
patch_variant_filters.py

Retune the variant browser's button-group filters for adaptive-sampling depth.

ALT count
---------
The thresholds were >10, >15, >20, which suit a deep hybrid-capture panel. On
this assay the on-target median runs closer to 10x, so a heterozygous variant
carries roughly 4-7 alt reads and every one of those buttons returns an empty
list. They become >1, >2, >5, which spans the range this data actually
occupies.

Callers
-------
Removed. The field is VariantCaller_Count and the options are >2, >3, >4.
Somatic SNV calling here uses a single caller, so the field is never populated
and no option can ever match. The sort entry that pairs with it goes too.

The filter is removed rather than repopulated with a constant: asserting a
caller count that was never measured would put a number in a clinical report
that no part of the pipeline produced.

Safety
------
Anchors are validated before anything is written, a timestamped .bak is kept,
and each edit is skipped if already applied. --dry-run writes nothing.

Usage
-----
    python3 bin/patch_variant_filters.py --dry-run
    python3 bin/patch_variant_filters.py
    python3 bin/patch_variant_filters.py --alt-thresholds 1,3,8
"""

import argparse
import os
import re
import shutil
import sys
from datetime import datetime


OLD_ALT_OPTIONS = re.compile(
    r'(\{ id: "gt10".*?\},\s*\n\s*\{ id: "gt15".*?\},\s*\n\s*\{ id: "gt20".*?\},)',
    re.DOTALL,
)


class PatchError(Exception):
    pass


def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def alt_options(thresholds, indent):
    lines = []
    for value in thresholds:
        label = ">%g" % value
        lines.append(
            '%s{ id: "gt%g", label: "%s",  test: function (v) '
            "{ return v !== null && v > %g; } },"
            % (indent, value, label, value)
        )
    return "\n".join(lines).lstrip()


def remove_object_block(text, key, description):
    """Remove `key: { ... },` including its trailing comma and newline."""
    index = text.find(key + ":")
    if index == -1:
        return text, False
    brace = text.find("{", index)
    if brace == -1:
        raise PatchError("no opening brace for %s" % description)
    depth = 0
    position = brace
    while position < len(text):
        char = text[position]
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                end = position + 1
                if text[end : end + 1] == ",":
                    end += 1
                if text[end : end + 1] == "\n":
                    end += 1
                line_start = text.rfind("\n", 0, index) + 1
                if not text[line_start:index].strip():
                    index = line_start
                return text[:index] + text[end:], True
        position += 1
    raise PatchError("unbalanced braces in %s" % description)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", default=".")
    parser.add_argument(
        "--alt-thresholds",
        default="1,2,5",
        help="comma-separated ALT count thresholds (default 1,2,5)",
    )
    parser.add_argument(
        "--keep-callers",
        action="store_true",
        help="leave the Callers filter in place",
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    path = os.path.join(
        args.repo, "bin", "dashboard_builder", "assets", "js", "variant-browser.js"
    )
    if not os.path.isfile(path):
        print("ERROR: not found: %s" % path, file=sys.stderr)
        return 2

    try:
        thresholds = [float(t) for t in args.alt_thresholds.split(",") if t.strip()]
    except ValueError:
        print("ERROR: --alt-thresholds must be numbers", file=sys.stderr)
        return 2
    if not thresholds:
        print("ERROR: at least one threshold is required", file=sys.stderr)
        return 2

    text = read(path)
    actions = []

    match = OLD_ALT_OPTIONS.search(text)
    if match:
        indent = " " * 8
        text = text[: match.start()] + alt_options(thresholds, indent) + text[match.end():]
        actions.append(
            "EDIT  ALT count thresholds -> %s"
            % ", ".join(">%g" % t for t in thresholds)
        )
    elif 'id: "gt%g"' % thresholds[0] in text:
        actions.append("SKIP  ALT count thresholds: already patched")
    else:
        print("ERROR: ALT count option block not found in the expected form.",
              file=sys.stderr)
        print("       Nothing was written.", file=sys.stderr)
        return 2

    if not args.keep_callers:
        try:
            text, removed = remove_object_block(
                text, "callers_buttons", "Callers filter"
            )
        except PatchError as error:
            print("ERROR: %s" % error, file=sys.stderr)
            return 2
        if removed:
            actions.append("EDIT  remove the Callers filter")
        else:
            actions.append("SKIP  Callers filter: already removed")

        # The sort entry that pairs with it.
        index = text.find('{ id: "callers_desc"')
        if index != -1:
            end = text.find("\n", index)
            line_start = text.rfind("\n", 0, index) + 1
            text = text[:line_start] + text[end + 1 :]
            actions.append("EDIT  remove the Callers sort option")

    for action in actions:
        print(action)

    if not [a for a in actions if a.startswith("EDIT")]:
        print("\nNothing to do.")
        return 0

    if args.dry_run:
        print("\nDry run: nothing written.")
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy2(path, "%s.bak_filters_%s" % (path, stamp))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    print("\nApplied. Backup written alongside.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
