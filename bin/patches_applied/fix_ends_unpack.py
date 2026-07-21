#!/usr/bin/env python3
"""
fix_ends_unpack.py
==================
Follow-up to apply_cytoband_partner_annotation.py.

The cytoband patch extended the _ends tuples in merge_translocations.py from
3 elements (chrom, pos, gene) to 4 (chrom, pos, gene, source). One consumer,
ig_aware_union(), unpacks each end into exactly three names:

    (ca, pa, ga), (cb, pb, gb) = rep["_ends"]

which raises "too many values to unpack (expected 3)" against a 4-tuple.
ig_aware_union clusters on position and anchor gene only; the source token is
not needed there, so we ignore the 4th element:

    (ca, pa, ga, _), (cb, pb, gb, _) = rep["_ends"]

Audit of all _ends consumers (verified against the file): only this one line
unpacks the inner tuple into a fixed arity. All other references either unpack
the *pair* of ends (end1, end2 = ...), or index (end[2]/end[3]), both of which
tolerate the 4-tuple.

Idempotent: writes a .bak, refuses to apply twice.

Usage (from repo root):
  python3 fix_ends_unpack.py --repo /goast/mm-awgs-nextflow
"""
import argparse
import datetime as _dt
import os
import shutil
import sys

OLD = '        (ca, pa, ga), (cb, pb, gb) = rep["_ends"]'
NEW = '        (ca, pa, ga, _), (cb, pb, gb, _) = rep["_ends"]'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True, help="Path to mm-awgs-nextflow repo root.")
    args = ap.parse_args()

    path = os.path.join(args.repo, "bin", "merge_translocations.py")
    if not os.path.isfile(path):
        sys.stderr.write(f"ERROR: not found: {path}\n")
        sys.exit(1)

    with open(path) as fh:
        text = fh.read()

    if NEW in text:
        sys.stderr.write(f"SKIP (already fixed): {path}\n")
        return

    if OLD not in text:
        sys.stderr.write(
            f"ERROR: expected line not found in {path}:\n  {OLD!r}\n"
            f"The file differs from the reviewed version; not modifying.\n")
        sys.exit(1)

    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.bak_endsfix_{stamp}"
    shutil.copy2(path, bak)

    text = text.replace(OLD, NEW, 1)
    with open(path, "w") as fh:
        fh.write(text)

    sys.stderr.write(f"PATCHED: {path}\n  backup: {bak}\n")
    sys.stderr.write("  Fixed ig_aware_union _ends unpack (3 -> 4 tuple).\n")


if __name__ == "__main__":
    main()
