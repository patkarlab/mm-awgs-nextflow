#!/usr/bin/env python3
"""
patch_cohort_index.py

Replace the hybrid-capture columns in the cohort index with metrics this assay
produces.

Removed
-------
Mean cov (mean of per-exon means), % >= 100x, Fold-80, % dup, FLT3-ITD.

All five come from Picard HsMetrics or a per-exon coverage table, neither of
which exists here, so every cell rendered blank.

Added
-----
Median depth, regions below 10x, rearrangements.

Median rather than mean is deliberate. The panel mixes focal gene-body windows
with megabase breakpoint windows, so a panel-wide mean is pulled down by the
wide ones and understates the focal regions. The median across regions is the
representative figure, and the low-coverage count sits beside it so a good
median cannot hide a set of regions with nothing in them.

Safety
------
Anchors are validated before anything is written, a timestamped .bak is kept,
and the edit is skipped if already applied. --dry-run writes nothing.

Usage
-----
    python3 bin/patch_cohort_index.py --dry-run
    python3 bin/patch_cohort_index.py
"""

import argparse
import os
import re
import shutil
import sys
from datetime import datetime


OLD_HEADERS = [
    "<th>Mean cov</th>",
    "<th>% &ge; 100x</th>",
    "<th>Fold-80</th>",
    "<th>% dup</th>",
    "<th>FLT3-ITD</th>",
]

NEW_HEADERS = """          <th>Median depth</th>
          <th>Regions &lt;10x</th>
          <th>Rearrangements</th>
"""

NEW_CELLS = """            <td>{{ "%.1f"|format(s.qc.depth.median) ~ "x" if s.qc and s.qc.depth else '' }}</td>
            <td>{{ s.qc.depth.n_below_10 if s.qc and s.qc.depth else '' }}</td>
            <td>{{ s.translocations.n_translocations if s.translocations and s.translocations.found else '' }}</td>
"""

SENTINEL = "s.qc.depth.median"


class PatchError(Exception):
    pass


def read(path):
    with open(path, encoding="utf-8") as handle:
        return handle.read()


def drop_line_containing(text, needle, description):
    """Remove the whole line holding needle."""
    index = text.find(needle)
    if index == -1:
        raise PatchError("anchor not found: %s (%r)" % (description, needle))
    start = text.rfind("\n", 0, index) + 1
    end = text.find("\n", index)
    end = len(text) if end == -1 else end + 1
    return text[:start] + text[end:]


def replace_between(text, after, before, replacement, description):
    """Replace everything between the end of `after` and the start of `before`."""
    start = text.find(after)
    if start == -1:
        raise PatchError("anchor not found: %s (start)" % description)
    start += len(after)
    end = text.find(before, start)
    if end == -1:
        raise PatchError("anchor not found: %s (end)" % description)
    line_start = text.rfind("\n", 0, end) + 1
    return text[:start] + "\n" + replacement + text[line_start:]


def drop_flt3_cell(text):
    """Remove the <td> holding the FLT3 conditional block."""
    index = text.find("{%- if s.flt3")
    if index == -1:
        return text, False
    start = text.rfind("<td>", 0, index)
    if start == -1:
        raise PatchError("no <td> before the FLT3 block")
    end = text.find("</td>", index)
    if end == -1:
        raise PatchError("no </td> after the FLT3 block")
    end += len("</td>")
    line_start = text.rfind("\n", 0, start) + 1
    if not text[line_start:start].strip():
        start = line_start
    if text[end : end + 1] == "\n":
        end += 1
    return text[:start] + text[end:], True


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--repo", default=".")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    path = os.path.join(
        args.repo, "bin", "dashboard_builder", "templates", "cohort_index.html.j2"
    )
    if not os.path.isfile(path):
        print("ERROR: not found: %s" % path, file=sys.stderr)
        return 2

    text = read(path)
    if SENTINEL in text:
        print("SKIP  cohort index: already patched")
        return 0

    actions = []
    try:
        # Header: first hybrid-capture column becomes the new set, the rest go.
        index = text.find(OLD_HEADERS[0])
        if index == -1:
            raise PatchError("anchor not found: Mean cov header")
        line_start = text.rfind("\n", 0, index) + 1
        line_end = text.find("\n", index) + 1
        text = text[:line_start] + NEW_HEADERS + text[line_end:]
        actions.append("EDIT  header: Mean cov -> median depth, low-coverage, rearrangements")

        for header in OLD_HEADERS[1:]:
            text = drop_line_containing(text, header, "header %s" % header)
            actions.append("EDIT  header: remove %s" % header)

        # Body: the metric cells sit between the sample link and the clinical
        # variant count, so both ends are anchored rather than the cells named.
        text = replace_between(
            text,
            "{{ s.sample }}</a></td>",
            "{{ s.clinical.n if s.clinical else '' }}",
            NEW_CELLS,
            "metric cells",
        )
        actions.append("EDIT  body: replace the four capture-metric cells")

        text, removed = drop_flt3_cell(text)
        if removed:
            actions.append("EDIT  body: remove the FLT3 cell")

    except PatchError as error:
        print("ERROR: %s" % error, file=sys.stderr)
        print("Nothing was written.", file=sys.stderr)
        return 2

    for action in actions:
        print(action)

    if args.dry_run:
        print("\nDry run: %d edit(s) planned, nothing written." % len(actions))
        return 0

    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    shutil.copy2(path, "%s.bak_cohort_%s" % (path, stamp))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(text)
    print("\nApplied %d edit(s). Backup written alongside." % len(actions))
    return 0


if __name__ == "__main__":
    sys.exit(main())
