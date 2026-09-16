#!/usr/bin/env python3
"""
collapse_bnd_mates.py

Collapse each pair of VCF breakend (BND) mates to a single record before the
SURVIVOR merge, so every caller contributes one record per junction.

Why
---
Sniffles and CuteSV write one BND record per junction. Severus (and SAVANA)
write both mates, linked by MATEID, and each mate enters SURVIVOR as its own
call. The two mates of one Severus junction can then land in two different
SURVIVOR clusters, each carrying a different SUPP_VEC, and the annotated
table shows two rows for one event with the caller set split between them.
merge_translocations.py repairs most of that positionally, downstream and
after the fact; this removes the cause.

Which mate is kept
------------------
The one whose primary position (CHROM, POS) sorts first in the reference
contig order taken from the VCF header, then by position. That is the
convention the single-record callers use, so a Severus junction enters the
merge in the same orientation as its Sniffles or CuteSV counterpart. The
mate is dropped; nothing about the kept record is rewritten.

What is not touched
-------------------
Records without MATEID, any non-BND record, and a BND whose mate is absent
from the file (its MATEID points at nothing) all pass through unchanged. The
script is therefore a no-op on Sniffles and CuteSV output and safe to apply
to every caller uniformly.

Standard library only. Reads plain or bgzipped VCF; writes plain VCF to
stdout or --output.

Usage:
  collapse_bnd_mates.py <in.vcf[.gz]> [--output out.vcf]
"""

import argparse
import gzip
import sys


def open_any(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def info_value(info, *keys):
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            if k in keys:
                return v
    return None


def fallback_chrom_key(chrom):
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(name.upper(), 99), name.upper())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("vcf")
    ap.add_argument("--output", default=None, help="write here instead of stdout")
    args = ap.parse_args()

    header = []
    records = []          # (line, chrom, pos, id, mateid, is_bnd)
    contig_order = {}
    with open_any(args.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                header.append(line)
                if line.startswith("##contig=<") and "ID=" in line:
                    cid = line.split("ID=", 1)[1].split(",", 1)[0].split(">", 1)[0]
                    contig_order.setdefault(cid, len(contig_order))
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            info = cols[7]
            svtype = info_value(info, "SVTYPE")
            mateid = info_value(info, "MATEID", "MATE_ID") if svtype == "BND" else None
            try:
                pos = int(cols[1])
            except ValueError:
                pos = 0
            records.append((line, cols[0], pos, cols[2], mateid, svtype == "BND"))

    def chrom_rank(chrom):
        if chrom in contig_order:
            return (0, contig_order[chrom], "")
        return (1,) + fallback_chrom_key(chrom)

    by_id = {r[3]: r for r in records}

    dropped = set()
    n_pairs = n_orphans = 0
    for line, chrom, pos, rid, mateid, is_bnd in records:
        if not is_bnd or not mateid or rid in dropped:
            continue
        # A MATEID may list several mates; take the first that exists.
        mate = None
        for mid in mateid.split(","):
            if mid in by_id and mid != rid:
                mate = by_id[mid]
                break
        if mate is None:
            n_orphans += 1
            continue
        if mate[3] in dropped:
            continue
        me_key = (chrom_rank(chrom), pos, rid)
        mate_key = (chrom_rank(mate[1]), mate[2], mate[3])
        loser = mate[3] if me_key <= mate_key else rid
        dropped.add(loser)
        n_pairs += 1

    out = open(args.output, "w") if args.output else sys.stdout
    try:
        for h in header:
            out.write(h)
        kept = 0
        for line, chrom, pos, rid, mateid, is_bnd in records:
            if rid in dropped:
                continue
            out.write(line)
            kept += 1
    finally:
        if args.output:
            out.close()

    sys.stderr.write(
        f"collapse_bnd_mates: {len(records)} records in, {kept} out; "
        f"{n_pairs} mate pairs collapsed, {n_orphans} BND with absent mate kept\n")


if __name__ == "__main__":
    main()
