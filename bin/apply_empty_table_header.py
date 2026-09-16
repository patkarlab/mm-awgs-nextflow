#!/usr/bin/env python3
"""
apply_empty_table_header.py

When filter_v6_somatic_candidates.py keeps no rows for a sample, its column
list is derived from the kept records and is therefore empty, so
csv.DictWriter writes a blank line and no header. alias_variant_table.py
passes that emptiness through as a zero-byte file, the dashboard gets a
table with no columns, and the strict bundle check (correctly) refuses it.
Seen on 11F20265533: 22 candidates, 18 off-panel, 4 in-panel and none
protein-altering, 0 kept.

Fix: when there are no kept rows, take the column order from PREFERRED_COLS
plus the input file's own columns, so an empty result is still a valid table
with a header. Nothing changes for non-empty results.

Idempotent; sentinel `empty_table_header_v1`; --check; --revert.
"""

import argparse
import os
import shutil
import sys

SENTINEL = "empty_table_header_v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = "bin/filter_v6_somatic_candidates.py"

OLD = (
    "    # Column order.\n"
    "    all_cols = set()\n"
    "    for r in kept:\n"
    "        all_cols.update(r.keys())\n"
    "    ordered = [c for c in PREFERRED_COLS if c in all_cols]\n"
    "    ordered += [c for c in all_cols if c not in ordered]\n"
)
NEW = (
    "    # Column order.\n"
    "    #\n"
    "    # empty_table_header_v1: with no kept rows the column set derived from\n"
    "    # the records is empty and DictWriter writes a blank line, not a header.\n"
    "    # Downstream then sees a zero-byte table. An empty result is a legitimate\n"
    "    # outcome and must still be a valid table, so the header falls back to\n"
    "    # the preferred columns plus whatever the input carried.\n"
    "    all_cols = set()\n"
    "    for r in kept:\n"
    "        all_cols.update(r.keys())\n"
    "    if not kept:\n"
    "        all_cols.update(PREFERRED_COLS)\n"
    "        all_cols.update(in_cols)\n"
    "    ordered = [c for c in PREFERRED_COLS if c in all_cols]\n"
    "    ordered += [c for c in all_cols if c not in ordered]\n"
)


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true")
    g.add_argument("--revert", action="store_true")
    args = ap.parse_args()

    path = os.path.join(REPO, TARGET)
    if args.revert:
        if os.path.isfile(path + ".bak"):
            shutil.move(path + ".bak", path)
            print(f"restored {TARGET}")
        else:
            print(f"no .bak  {TARGET}")
        return

    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    if SENTINEL in text:
        print(f"APPLIED  {TARGET}")
        return
    if text.count(OLD) != 1:
        print(f"NO-MATCH {TARGET}: anchor found {text.count(OLD)} times")
        sys.exit(1)
    if args.check:
        print(f"READY    {TARGET}")
        return
    shutil.copy2(path, path + ".bak")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text.replace(OLD, NEW, 1))
    print(f"patched  {TARGET} (backup: {TARGET}.bak)")


if __name__ == "__main__":
    main()
