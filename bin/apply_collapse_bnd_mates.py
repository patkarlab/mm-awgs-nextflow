#!/usr/bin/env python3
"""
apply_collapse_bnd_mates.py

Handoff open item 7 (2026-09-15): collapse SV mates on MATEID.

SURVIVOR rewrites every breakend as SVTYPE=TRA and drops MATEID, so there is
nothing to collapse on after the merge; the collapse has to happen before.
Severus writes both mates of every breakend (198 records = 99 junctions on
11F20265231, every one with MATE_ID); Sniffles and CuteSV write one record
per junction. Each Severus junction therefore enters the merge twice, the
two mates can land in different SURVIVOR clusters with different SUPP_VEC,
and merge_translocations.py repairs it positionally, after the fact.

Change: SURVIVOR_MERGE runs bin/collapse_bnd_mates.py on all four caller
VCFs before writing vcflist.txt. The script keeps the mate whose primary
position is canonical-first (header contig order, then position) — the
convention Sniffles and CuteSV already use — and is a no-op on a file with
no MATEIDs. Published caller VCFs are untouched; augment_sv_support.py reads
those and matches orientation-agnostically.

Risk, and the check that covers it: if "canonical-first" is not the
orientation Sniffles/CuteSV use, Severus support would stop merging with
theirs and 'severus' would vanish from the callers column on junctions it
currently supports. bin/check_caller_retention.sh compares, per sample, the
number of merged junctions carrying each caller before and after. Run it
before committing.

Idempotent; sentinel `collapse_bnd_mates_v1`; --check; --revert.
"""

import argparse
import os
import shutil
import sys

SENTINEL = "collapse_bnd_mates_v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = "bin/collapse_bnd_mates.py"

EDITS = [("modules/local/survivor_merge.nf", [(
    "    # The caller order written to vcflist.txt is the order SURVIVOR uses for\n",
    "    # collapse_bnd_mates_v1. Severus and SAVANA write both mates of each\n"
    "    # breakend; Sniffles and CuteSV write one record per junction. Collapse\n"
    "    # every caller to one record per junction, keeping the mate on the\n"
    "    # canonical-first end, so a junction enters the merge once per caller\n"
    "    # and SUPP_VEC cannot be split across its two mates. No-op on a file\n"
    "    # without MATEIDs. The published caller VCFs are not changed.\n"
    "    for c in sniffles cutesv severus savana; do\n"
    "        collapse_bnd_mates.py \\${c}.vcf --output \\${c}.collapsed.vcf\n"
    "        mv \\${c}.collapsed.vcf \\${c}.vcf\n"
    "    done\n"
    "\n"
    "    # The caller order written to vcflist.txt is the order SURVIVOR uses for\n",
)])]


def read(path):
    with open(path, encoding="utf-8") as fh:
        return fh.read()


def check():
    problems = 0
    if os.path.isfile(os.path.join(REPO, SCRIPT)):
        print(f"PRESENT  {SCRIPT}")
    else:
        print(f"MISSING  {SCRIPT} (copy it from the inbox first)"); problems += 1
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        text = read(path)
        if SENTINEL in text:
            print(f"APPLIED  {rel}"); continue
        bad = [old for old, _ in pairs if text.count(old) != 1]
        if bad:
            print(f"NO-MATCH {rel}"); problems += 1
        else:
            print(f"READY    {rel}")
    return 1 if problems else 0


def apply():
    if not os.path.isfile(os.path.join(REPO, SCRIPT)):
        sys.exit(f"ERROR: {SCRIPT} not found; nothing changed.")
    for rel, pairs in EDITS:
        path = os.path.join(REPO, rel)
        text = read(path)
        if SENTINEL in text:
            print(f"skip     {rel} (already applied)"); continue
        for old, new in pairs:
            if text.count(old) != 1:
                sys.exit(f"ERROR: anchor not found exactly once in {rel}; nothing changed.")
            text = text.replace(old, new, 1)
        shutil.copy2(path, path + ".bak")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(text)
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
