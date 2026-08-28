#!/usr/bin/env python3
"""
derive_excluded_junctions.py
============================

Propose artefact junctions from a cohort's own annotated tables.

The criterion is coordinate identity across unrelated patients, not
recurrence. Somatic breakpoints do not recur to the nucleotide between
individuals: repair at a real junction is imprecise, so two patients sharing
a rearrangement share the intron, not the base. Agreement to within a few
bases across unrelated samples is a property of the reference or the aligner,
not of the tumours.

Recurrence alone must never become the criterion, and this script cannot be
used to make it one. Plasma cell neoplasms are defined by recurrent events:
t(11;14) is expected in 15-20% of patients and t(4;14) in 10-15%. A list
built on "seen in several samples" would exclude the panel's entire purpose.
Two guards enforce the distinction:

  A junction whose row carries a dictionary name or a tier is never
  proposed, whatever its coordinates do. That mirrors the override in the
  annotator, which is not configurable either.

  The default tolerance is 200 bp on BOTH breakends. A real shared
  rearrangement will not satisfy that across unrelated patients; a mapping
  artefact will.

Output is a proposal, not a decision. Every row is written commented out,
with the samples and coordinates that produced it, for a human to read and
uncomment. Nothing is excluded until someone does.

Usage
-----
  derive_excluded_junctions.py \\
      --annotated results/t2t/calls/mm_annotated/*.mm_annotated.tsv \\
      --min-samples 3 \\
      --tolerance 200 \\
      --output assets/mm_excluded_junctions.proposed.tsv

--min-samples is the number of distinct samples that must carry the junction
at coordinate identity. Three is a deliberate floor: two samples agreeing can
happen, and a third makes chance implausible.

Timepoints of one sample are not independent. Pass --sample-groups to declare
which identifiers belong to the same patient, one group per line, tab or
comma separated; grouped identifiers are counted once.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

__version__ = "0.1.0"


def load_groups(path):
    """Sample id -> group id. Ungrouped ids are their own group."""
    groups = {}
    if path is None:
        return groups
    with open(path) as fh:
        for i, line in enumerate(fh):
            if not line.strip() or line.startswith("#"):
                continue
            ids = [t.strip() for t in line.replace(",", "\t").split("\t")
                   if t.strip()]
            for sid in ids:
                groups[sid] = f"group{i}"
    return groups


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotated", required=True, nargs="+",
                    help="One or more *.mm_annotated.tsv files.")
    ap.add_argument("--min-samples", type=int, default=3)
    ap.add_argument("--tolerance", type=int, default=200,
                    help="Bases either breakend may differ and still count as "
                         "the same junction [200].")
    ap.add_argument("--sample-groups", default=None, type=Path,
                    help="Identifiers belonging to one patient, one group per "
                         "line. Grouped identifiers count once.")
    ap.add_argument("--include-intrachromosomal", action="store_true",
                    help="Also consider same-chromosome events. Off by "
                         "default: they dominate the callset and recurrent "
                         "small indels are a population-frequency problem "
                         "rather than one to enumerate by coordinate.")
    ap.add_argument("--partner-tol", type=int, default=None,
                    help="Tolerance on the side that clusters tightly, when "
                         "the other side is dispersed. Defaults to "
                         "--tolerance. An Ig breakend scatters over the "
                         "locus while its partner does not, so requiring "
                         "both sides within one tolerance splits a single "
                         "junction into several proposals.")
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    groups = load_groups(args.sample_groups)
    include_intra = args.include_intrachromosomal

    # Canonical orientation, so a reciprocal mate lands in the same bucket.
    obs = defaultdict(list)
    protected = 0
    n_rows = 0
    for path in args.annotated:
        with open(path) as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n_rows += 1
                if (r.get("known_mm_pair") or "").strip() or \
                   (r.get("tier") or "").strip():
                    protected += 1
                    continue
                ca, cb = r.get("chrom_a"), r.get("chrom_b")
                # Inter-chromosomal only, by default. Same-chromosome DEL
                # and INS dominate the callset: on one cohort they were
                # 8868 of 9735 rows and produced 1019 proposals where 42
                # were junctions. Recurrent small indels at shared
                # coordinates are germline variation and alignment noise,
                # which is a population-frequency problem, not something to
                # enumerate by coordinate.
                if not include_intra and ca == cb:
                    continue
                try:
                    pa, pb = int(r.get("pos_a")), int(r.get("pos_b"))
                except (TypeError, ValueError):
                    continue
                if not ca or not cb:
                    continue
                if (ca, pa) > (cb, pb):
                    ca, pa, cb, pb = cb, pb, ca, pa
                sid = (r.get("sample") or Path(path).name).strip()
                obs[(ca, cb)].append((pa, pb, groups.get(sid, sid), sid,
                                      r.get("gene_a", ""), r.get("gene_b", "")))

    ptol = args.partner_tol if args.partner_tol is not None else args.tolerance
    proposals = []
    for (ca, cb), pts in obs.items():
        pts.sort()
        used = [False] * len(pts)
        for i, (pa, pb, g, sid, ga, gb) in enumerate(pts):
            if used[i]:
                continue
            cluster = [(pa, pb, g, sid, ga, gb)]
            used[i] = True
            for j in range(i + 1, len(pts)):
                if used[j]:
                    continue
                qa, qb, g2, sid2, ga2, gb2 = pts[j]
                # The scan is sorted on side a, so it can stop once side
                # a is out of range -- but the range is the wider of the
                # two tolerances, since a match may allow side a to be
                # the dispersed one.
                if qa - pa > max(args.tolerance, ptol):
                    break
                # Asymmetric. One side may scatter across a locus while
                # the other is tight: six IGH proposals in one run were a
                # single junction whose partner clustered within 1.3 kb
                # while the IGH side spread over 640 kb. A match needs one
                # side within --partner-tol and the other within
                # --tolerance, either way round.
                da, db = abs(qa - pa), abs(qb - pb)
                if ((da <= ptol and db <= args.tolerance) or
                        (da <= args.tolerance and db <= ptol)):
                    cluster.append(pts[j])
                    used[j] = True
            indep = {c[2] for c in cluster}
            if len(indep) < args.min_samples:
                continue
            sas = sorted(c[0] for c in cluster)
            sbs = sorted(c[1] for c in cluster)
            proposals.append({
                "chrom_a": ca, "pos_a": sas[len(sas) // 2],
                "chrom_b": cb, "pos_b": sbs[len(sbs) // 2],
                "tolerance": args.tolerance,
                "n_indep": len(indep),
                # A count, not identifiers. Sequencing IDs are lab case numbers
                # and this file is committed to a public repository. The
                # identifiers remain in the run's stderr, which is not.
                "seen_in": str(len({c[3] for c in cluster})),
                "genes": f"{cluster[0][4]} x {cluster[0][5]}",
                "spread": f"{max(sas) - min(sas)}/{max(sbs) - min(sbs)} bp",
            })

    proposals.sort(key=lambda p: (-p["n_indep"], p["chrom_a"], p["pos_a"]))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w") as out:
        out.write(
            "# Proposed artefact junctions, derived from this cohort.\n"
            "#\n"
            "# EVERY ROW IS COMMENTED OUT. Nothing here is excluded until a\n"
            "# human reads it and uncomments it. Read each one against the\n"
            "# alignment before you do.\n"
            "#\n"
            f"# Criterion: coordinate identity within {args.tolerance} bp on\n"
            f"# BOTH breakends, in at least {args.min_samples} independent\n"
            "# samples. Somatic breakpoints do not recur to the nucleotide\n"
            "# between individuals, so agreement this tight is a property of\n"
            "# the reference or the aligner rather than of the tumours.\n"
            "#\n"
            "# Recurrence alone is NOT the criterion and must not become one.\n"
            "# This disease is defined by recurrent events: t(11;14) is\n"
            "# expected in 15-20% of patients. Rows carrying a dictionary\n"
            "# name or a tier were excluded from consideration before\n"
            "# clustering and cannot appear below.\n"
            "#\n"
            "# A junction on an Ig anchor deserves particular scrutiny: the\n"
            "# anchor rule makes any IGH, IGK, IGL or MYC junction reportable,\n"
            "# so an artefact there reaches the report on every sample.\n"
            "#\n"
            f"# Generated from {len(args.annotated)} table(s), {n_rows} rows,\n"
            f"# of which {protected} were protected as named or graded.\n"
            "#\n"
            "# CHROM_A\tPOS_A\tCHROM_B\tPOS_B\tTOLERANCE\tN_SAMPLES\tNOTE\n")
        for p in proposals:
            out.write(
                f"# {p['n_indep']} independent sample(s); breakpoint spread "
                f"{p['spread']}; {p['genes']}\n"
                f"#{p['chrom_a']}\t{p['pos_a']}\t{p['chrom_b']}\t{p['pos_b']}"
                f"\t{p['tolerance']}\t{p['seen_in']}\t\n")

    sys.stderr.write(
        f"{n_rows} row(s) across {len(args.annotated)} table(s); "
        f"{protected} protected as named or graded; "
        f"{len(proposals)} proposal(s) -> {args.output}\n")
    if not proposals:
        sys.stderr.write(
            "No junction met the criterion. That is a normal result and a\n"
            "good one: it means no coordinate-identical junction recurs\n"
            "across this cohort at the tolerance given.\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
