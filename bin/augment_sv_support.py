#!/usr/bin/env python3
"""
augment_sv_support.py

The SURVIVOR-merged VCF drops per-caller read support (it keeps only SUPP/
SUPP_VEC and coordinates), so support_reads cannot be recovered from the merged
file or from mm_annotated.tsv. This script layers the real per-caller support
back on, read directly from each caller's own VCF:

  Sniffles : INFO/SUPPORT
  CuteSV   : INFO/RE
  Severus  : FORMAT/DV   (fallback INFO/SUPP_READS first colon-field)

For every row in an mm_annotated.tsv it matches the junction (orientation-
agnostic, within --tol bp on BOTH breakpoints to absorb the few-bp shifts
SURVIVOR introduces) against each caller's calls and writes:

  support_sniffles, support_cutesv, support_severus   per-caller read counts
  support_reads (overwritten)                          = max of those present

Matching is by coordinate only; no specific variants or findings are hardcoded.
Same-chromosome SV types (DEL/DUP/INV/INS) are matched on the primary breakpoint
plus END; translocations (TRA/BND) on both breakends.

Standard-library only.

Usage:
  python3 augment_sv_support.py \\
      --annotated  <sample>.mm_annotated.tsv \\
      --sniffles   <sample>.sniffles.t2t.vcf.gz \\
      --cutesv     <sample>.cutesv.t2t.vcf.gz \\
      --severus    <sample>.severus.vcf \\
      --output     <sample>.mm_annotated.tsv      (in place is fine)
      [--tol 25]
"""

import argparse
import csv
import gzip
import os
import sys
from typing import Dict, List, Optional, Tuple


def open_any(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def chrom_key(chrom):
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(name.upper(), 99), name.upper())


def parse_info(info_field):
    d = {}
    for part in info_field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        elif part:
            d[part] = "True"
    return d


def parse_bnd_alt(alt):
    for opener, closer in (("]", "]"), ("[", "[")):
        if opener in alt:
            try:
                inner = alt.split(opener)[1].split(closer)[0]
                c, p = inner.rsplit(":", 1)
                return c, int(p)
            except (IndexError, ValueError):
                continue
    return None, None


def canon(ca, pa, cb, pb):
    """Canonical orientation-agnostic ordering of the two ends."""
    ka, kb = (chrom_key(ca), pa), (chrom_key(cb), pb)
    return (ca, pa, cb, pb) if ka <= kb else (cb, pb, ca, pa)


def caller_support(value_kind, info, fmt, sample):
    """Extract this caller's own support count from one of its records."""
    if value_kind == "RE":
        v = info.get("RE", "")
        return int(v) if str(v).isdigit() else None
    if value_kind == "SUPPORT":
        v = info.get("SUPPORT", "")
        return int(v) if str(v).isdigit() else None
    if value_kind == "TUMOUR_READ_SUPPORT":
        # SAVANA. An INFO field, like RE and SUPPORT, so no FORMAT lookup.
        v = info.get("TUMOUR_READ_SUPPORT", "")
        return int(v) if str(v).isdigit() else None
    if value_kind == "DV":
        if fmt and sample:
            f = dict(zip(fmt.split(":"), sample.split(":")))
            dv = f.get("DV", "")
            if dv.isdigit():
                return int(dv)
        sr = info.get("SUPP_READS", "")
        first = sr.split(":")[0] if sr else ""
        return int(first) if first.isdigit() else None
    return None


def load_caller(path, value_kind):
    """Return a list of (canon_a_chrom, a_pos, b_chrom, b_pos, support) for
    every SV record in a caller VCF."""
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open_any(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, _id, _ref, alt, _q, filt, info_field = cols[:8]
            pos = int(pos)
            info = parse_info(info_field)
            svtype = info.get("SVTYPE", "")
            fmt = cols[8] if len(cols) > 8 else ""
            sample = cols[9] if len(cols) > 9 else ""
            sup = caller_support(value_kind, info, fmt, sample)
            if sup is None:
                continue
            if svtype in ("BND", "TRA"):
                mc, mp = parse_bnd_alt(alt)
                if mc is None:
                    mc, end = info.get("CHR2"), info.get("END")
                    mp = int(end) if (end and end.isdigit()) else None
                if mc is None or mp is None:
                    continue
                ca, pa, cb, pb = canon(chrom, pos, mc, mp)
            else:
                end = info.get("END")
                if not (end and end.isdigit()):
                    continue
                ca, pa, cb, pb = canon(chrom, pos, chrom, int(end))
            out.append((ca, pa, cb, pb, sup, filt))
    return out


def load_nanomonsv(path):
    """nanomonsv_confirmation_v1: (canon ends, supporting reads, Is_Filter)
    for every row of a nanomonsv result.txt. Columns are read by header
    name; Dir_1/Dir_2 are not used because the ensemble match is
    orientation-agnostic. Filtered rows are kept so a junction nanomonsv
    saw and rejected can be told apart from one it never reported."""
    out = []
    if not path or not os.path.isfile(path):
        return out
    with open_any(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            try:
                c1, p1 = row["Chr_1"], int(row["Pos_1"])
                c2, p2 = row["Chr_2"], int(row["Pos_2"])
            except (KeyError, TypeError, ValueError):
                continue
            sup = row.get("Supporting_Read_Num_Tumor", "")
            sup = int(sup) if str(sup).isdigit() else None
            if sup is None:
                continue
            filt = (row.get("Is_Filter") or "").strip() or "PASS"
            ca, pa, cb, pb = canon(c1, p1, c2, p2)
            out.append((ca, pa, cb, pb, sup, filt))
    return out


def best_match(ca, pa, cb, pb, records, tol):
    """Support of the NEAREST caller record whose BOTH ends fall within tol bp
    (sum of per-end distances minimised). Nearest, not max, so that in dense
    loci where a caller emits several breakends a few bp apart, the count is
    attributed to the corresponding breakpoint rather than the strongest nearby
    one. Ties in distance are broken by higher support."""
    best_sup = None
    best_filt = None
    best_dist = None
    for (rca, rpa, rcb, rpb, sup, filt) in records:
        if rca == ca and rcb == cb and abs(rpa - pa) <= tol and abs(rpb - pb) <= tol:
            dist = abs(rpa - pa) + abs(rpb - pb)
            if (best_dist is None or dist < best_dist
                    or (dist == best_dist and (best_sup is None or sup > best_sup))):
                best_dist, best_sup, best_filt = dist, sup, filt
    return best_sup, best_filt


def main():
    ap = argparse.ArgumentParser(description="Layer per-caller read support onto mm_annotated.tsv.")
    ap.add_argument("--annotated", required=True)
    ap.add_argument("--sniffles", default=None)
    ap.add_argument("--cutesv", default=None)
    ap.add_argument("--severus", default=None)
    ap.add_argument("--savana", default=None,
                    help="SAVANA classified VCF (not classified.somatic.vcf; "
                         "the somatic model rejects low-support translocations "
                         "on this data)")
    ap.add_argument("--output", required=True)
    ap.add_argument("--tol", type=int, default=25, help="bp tolerance per breakpoint (default 25).")
    # nanomonsv_confirmation_v1
    ap.add_argument("--nanomonsv", default=None,
                    help="nanomonsv result.txt; adds support_nanomonsv, "
                         "filter_nanomonsv and nanomonsv_confirmed. Read beside "
                         "the ensemble, never folded into support_reads.")
    ap.add_argument("--nanomonsv-tol", type=int, default=100,
                    help="bp tolerance per breakend for the confirmation match "
                         "(default 100; SURVIVOR shifts positions, nanomonsv does not).")
    args = ap.parse_args()

    snf = load_caller(args.sniffles, "SUPPORT")
    cut = load_caller(args.cutesv, "RE")
    sev = load_caller(args.severus, "DV")
    sav = load_caller(args.savana, "TUMOUR_READ_SUPPORT")
    nano = load_nanomonsv(args.nanomonsv)
    nano_run = bool(args.nanomonsv)
    sys.stderr.write(
        f"loaded support records: sniffles={len(snf)} cutesv={len(cut)} "
        f"severus={len(sev)} savana={len(sav)} nanomonsv={len(nano)}"
        f"{'' if nano_run else ' (not run)'}\n")

    with open(args.annotated, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    # filter_savana carries SAVANA's somatic verdict without a column of its
    # own: SAVANA writes PASS or PREDICTED_NOISE into the VCF FILTER field.
    new_cols = ["support_sniffles", "support_cutesv", "support_severus",
                "support_savana",
                "filter_sniffles", "filter_cutesv", "filter_severus",
                "filter_savana",
                "filter_worst",
                "support_nanomonsv", "filter_nanomonsv", "nanomonsv_confirmed"]
    out_cols = in_cols + [c for c in new_cols if c not in in_cols]

    n_pop = 0
    for r in rows:
        ca_, pa_ = r.get("chrom_a", ""), r.get("pos_a", "")
        cb_, pb_ = r.get("chrom_b", ""), r.get("pos_b", "")
        if not (ca_ and pa_ and cb_ and pb_ and pa_.lstrip("-").isdigit() and pb_.lstrip("-").isdigit()):
            for c in new_cols:
                r[c] = ""
            continue
        ca, pa, cb, pb = canon(ca_, int(pa_), cb_, int(pb_))
        # Gate on the callers the merge already established for this junction
        # (from SUPP_VEC). Coordinate-matching only fetches that caller's read
        # count; it never credits a caller absent from this row, which prevents
        # a different nearby breakend in a dense locus from bleeding across.
        listed = {x.strip().lower() for x in (r.get("callers") or "").split(",") if x.strip()}
        def gated(name, recs):
            if listed and name not in listed:
                return (None, None)
            return best_match(ca, pa, cb, pb, recs, args.tol)
        s, sf = gated("sniffles", snf)
        c, cf = gated("cutesv", cut)
        v, vf = gated("severus", sev)
        a, af = gated("savana", sav)
        r["support_sniffles"] = str(s) if s is not None else ""
        r["support_cutesv"]   = str(c) if c is not None else ""
        r["support_severus"]  = str(v) if v is not None else ""
        r["support_savana"]   = str(a) if a is not None else ""
        r["filter_sniffles"]  = sf or ""
        r["filter_cutesv"]    = cf or ""
        r["filter_severus"]   = vf or ""
        r["filter_savana"]    = af or ""
        present = [x for x in (s, c, v, a) if x is not None]
        if present:
            r["support_reads"] = str(max(present))
            n_pop += 1
        # Worst FILTER across the callers that saw this junction.
        # SURVIVOR stamps PASS on every merged record, so the merged
        # filter column carries no quality information; this one does.
        # It is a stratifier, not a filter: a FISH-confirmed
        # rearrangement can carry a caller's low-support verdict and
        # still be real, so nothing downstream drops on it.
        seen = [f for f in (sf, cf, vf) if f]
        nonpass = sorted({f for f in seen if f != "PASS"})
        r["filter_worst"] = ";".join(nonpass) if nonpass else ("PASS" if seen else "")
        # nanomonsv_confirmation_v1. Not gated on the callers column: the
        # point is an independent method seeing the same junction. 'yes'
        # needs a PASS nanomonsv junction within tolerance on both ends;
        # a filtered nanomonsv hit is 'no' with its filter recorded; blank
        # means nanomonsv did not run for this sample.
        if nano_run:
            ns, nf = best_match(ca, pa, cb, pb, nano, args.nanomonsv_tol)
            r["support_nanomonsv"] = str(ns) if ns is not None else ""
            r["filter_nanomonsv"] = nf or ""
            r["nanomonsv_confirmed"] = "yes" if (ns is not None and nf == "PASS") else "no"
        else:
            r["support_nanomonsv"] = r["filter_nanomonsv"] = r["nanomonsv_confirmed"] = ""

    with open(args.output, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore", restval="", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)

    sys.stderr.write(f"populated support_reads on {n_pop}/{len(rows)} rows -> {args.output}\n")


if __name__ == "__main__":
    main()
