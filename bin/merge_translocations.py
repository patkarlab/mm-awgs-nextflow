#!/usr/bin/env python3
"""
merge_translocations.py

Collapse duplicate translocation (TRA) calls in an mm_annotated.tsv where the
same junction is reported multiple times - as reciprocal mates (Severus _1/_2
with A/B swapped) or as near-identical calls from different callers a few bp
apart. Two TRA rows are united when, after orientation normalisation, they share
the same chromosome pair AND BOTH breakpoints lie within --max-dist bp
(default 100). Non-TRA rows pass through unchanged.

Standard-library only; runs in any python3 (incl. awgs_sv).

Merge semantics per cluster:
  - orientation normalised: each row's two ends are ordered canonically
    (chrom-sorted), so a mate's swapped A/B collapses onto its partner;
  - representative row = most callers, then most support_reads, then first;
  - positions/genes taken from the representative (gene stays with its end;
    a panel-gene annotation is preferred over OFF_PANEL);
  - callers = union across members; n_callers recomputed from that union;
  - known_mm_pair / known_freq = first non-empty across members;
  - filter = PASS if any member is PASS;
  - two columns added: n_merged, merged_sv_ids.

Single-linkage clustering (transitive) is used; at 100 bp clusters are tiny so
chaining is not a concern.

Usage:
  python3 merge_translocations.py --input <sample>.mm_annotated.tsv [...] \
      [--outdir DIR] [--max-dist 100]
"""

import argparse
import csv
import os
import sys

CALLER_ORDER = {"Sniffles": 0, "CuteSV": 1, "Severus": 2, "nanomonsv": 3}


def chrom_sort_key(chrom):
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    return (1, {"X": 0, "Y": 1, "M": 2, "MT": 2}.get(name.upper(), 99), name.upper())


def to_int(value, default=None):
    try:
        return int(str(value).strip())
    except (ValueError, TypeError):
        return default


def canonical_ends(row):
    """Return ((chromA,posA,geneA),(chromB,posB,geneB)) chrom-sorted so that
    reciprocal mates normalise identically."""
    a = (row["chrom_a"], to_int(row["pos_a"], 0), row.get("gene_a", ""))
    b = (row["chrom_b"], to_int(row["pos_b"], 0), row.get("gene_b", ""))
    ka = (chrom_sort_key(a[0]), a[1])
    kb = (chrom_sort_key(b[0]), b[1])
    return (a, b) if ka <= kb else (b, a)


def callers_set(row):
    raw = (row.get("callers") or "").strip()
    return {c.strip() for c in raw.split(",") if c.strip()}


def best_gene(members, end_index):
    """Pick a gene for end 0 or 1: prefer a non-empty, non-OFF_PANEL value."""
    fallback = ""
    for m in members:
        g = m["_ends"][end_index][2].strip()
        if g and g != "OFF_PANEL":
            return g
        if g and not fallback:
            fallback = g
    return fallback


def first_nonempty(members, col):
    for m in members:
        v = (m.get(col) or "").strip()
        if v:
            return v
    return ""


def merge_cluster(members):
    # representative: most callers, then most support_reads, then first.
    def rep_key(m):
        return (to_int(m.get("n_callers"), 0),
                to_int(m.get("support_reads"), 0))
    rep = max(members, key=rep_key)

    end1, end2 = rep["_ends"]
    callers = set()
    for m in members:
        callers |= callers_set(m)
    callers_sorted = sorted(callers, key=lambda c: (CALLER_ORDER.get(c, 99), c))

    # max, not sum: the same junction seen by multiple callers is overlapping
    # read sets, so summing would double/triple-count the same molecules.
    supports = [to_int(m.get("support_reads")) for m in members]
    supports = [s for s in supports if s is not None]
    support_reads = str(max(supports)) if supports else ""

    def _max_col(col):
        vals = [to_int(m.get(col)) for m in members]
        vals = [v for v in vals if v is not None]
        return str(max(vals)) if vals else ""
    sup_snf = _max_col("support_sniffles")
    sup_cut = _max_col("support_cutesv")
    sup_sev = _max_col("support_severus")

    filt = "PASS" if any((m.get("filter") or "").strip() == "PASS" for m in members) \
        else (rep.get("filter") or "")

    return {
        "sample": rep.get("sample", ""),
        "sv_id": rep.get("sv_id", ""),
        "sv_type": "TRA",
        "filter": filt,
        "chrom_a": end1[0], "pos_a": str(end1[1]), "gene_a": best_gene(members, 0),
        "chrom_b": end2[0], "pos_b": str(end2[1]), "gene_b": best_gene(members, 1),
        "known_mm_pair": first_nonempty(members, "known_mm_pair"),
        "known_freq": first_nonempty(members, "known_freq"),
        "callers": ",".join(callers_sorted),
        "n_callers": str(len(callers_sorted)),
        "supp_vec": rep.get("supp_vec", ""),
        "support_reads": support_reads,
        "support_sniffles": sup_snf,
        "support_cutesv": sup_cut,
        "support_severus": sup_sev,
        "n_merged": str(len(members)),
        "merged_sv_ids": ",".join(m.get("sv_id", "") for m in members),
    }


def cluster_tra(tra_rows, max_dist):
    """Single-linkage union-find on canonical-pair + both-ends-within-threshold."""
    n = len(tra_rows)
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        parent[find(x)] = find(y)

    for i in range(n):
        e1i, e2i = tra_rows[i]["_ends"]
        for j in range(i + 1, n):
            e1j, e2j = tra_rows[j]["_ends"]
            if e1i[0] != e1j[0] or e2i[0] != e2j[0]:
                continue
            if abs(e1i[1] - e1j[1]) <= max_dist and abs(e2i[1] - e2j[1]) <= max_dist:
                union(i, j)

    clusters = {}
    for i in range(n):
        clusters.setdefault(find(i), []).append(tra_rows[i])
    # return clusters keyed by the smallest original index (preserve order)
    return clusters


def process_file(path, outdir, max_dist, keep_all_sv):
    stem = os.path.basename(path)
    # default: translocations-only file; --keep-all-sv: full table w/ TRA merged
    suffix = ".sv_tra_merged.tsv" if keep_all_sv else ".translocations.tsv"
    base = stem
    for s_ in (".mm_annotated.tsv", "_mm_annotated.tsv", ".tsv"):
        if base.endswith(s_):
            base = base[: -len(s_)]
            break
    out_name = base + suffix
    os.makedirs(outdir, exist_ok=True)
    out_path = os.path.join(outdir, out_name)

    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]

    out_cols = in_cols + [c for c in ("n_merged", "merged_sv_ids") if c not in in_cols]

    # Index TRA rows; precompute canonical ends.
    tra_idx = []
    for i, r in enumerate(rows):
        if (r.get("sv_type") or "").strip() == "TRA":
            r["_ends"] = canonical_ends(r)
            tra_idx.append(i)

    tra_rows = [rows[i] for i in tra_idx]
    clusters = cluster_tra(tra_rows, max_dist)

    # Map each TRA row's original position -> its cluster's representative output,
    # emitted once at the earliest member position.
    # Build cluster member lists keyed by min original index.
    pos_of = {id(r): tra_idx[k] for k, r in enumerate(tra_rows)}
    merged_by_first = {}
    member_first = {}
    for _, members in clusters.items():
        first_pos = min(pos_of[id(m)] for m in members)
        merged_by_first[first_pos] = merge_cluster(members)
        for m in members:
            member_first[pos_of[id(m)]] = first_pos

    out_rows = []
    n_tra_in = len(tra_rows)
    for i, r in enumerate(rows):
        if (r.get("sv_type") or "").strip() == "TRA":
            if member_first.get(i) == i:        # representative emission point
                out_rows.append(merged_by_first[i])
            # else: a non-first cluster member -> skip
        elif keep_all_sv:
            rr = {k: r.get(k, "") for k in in_cols}
            rr["n_merged"] = "1"
            rr["merged_sv_ids"] = r.get("sv_id", "")
            out_rows.append(rr)
        # else: non-TRA rows are dropped (translocations-only output)

    with open(out_path, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=out_cols, delimiter="\t",
                           extrasaction="ignore", restval="", lineterminator="\n")
        w.writeheader()
        for rr in out_rows:
            w.writerow(rr)

    n_tra_out = sum(1 for rr in out_rows if rr.get("sv_type") == "TRA")
    scope = "full SV table" if keep_all_sv else "translocations only"
    sys.stderr.write(
        f"[{stem}] TRA {n_tra_in} -> {n_tra_out} merged; "
        f"{len(out_rows)} rows written ({scope}); {out_path}\n"
    )
    return out_path


def main():
    ap = argparse.ArgumentParser(description="Unite near-identical translocation calls in mm_annotated.tsv.")
    ap.add_argument("-i", "--input", required=True, nargs="+", help="One or more *.mm_annotated.tsv files.")
    ap.add_argument("-o", "--outdir", default=None, help="Output dir (default: alongside each input).")
    ap.add_argument("--max-dist", type=int, default=100, help="Max bp between breakpoints to unite (default 100).")
    ap.add_argument("--keep-all-sv", action="store_true", help="Keep all SV types in output (TRA merged, others passthrough). Default: translocations only.")
    args = ap.parse_args()
    for path in args.input:
        if not os.path.isfile(path):
            sys.stderr.write(f"ERROR: not found: {path}\n"); sys.exit(1)
        outdir = args.outdir or os.path.dirname(os.path.abspath(path))
        process_file(path, outdir, args.max_dist, args.keep_all_sv)


if __name__ == "__main__":
    main()
