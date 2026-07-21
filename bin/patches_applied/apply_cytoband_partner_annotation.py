#!/usr/bin/env python3
"""
apply_cytoband_partner_annotation.py
====================================
Idempotent patch: characterize off-panel translocation partners by cytoband
instead of the bare "OFF_PANEL" token, across two scripts:

  bin/annotate_mm_translocations.py
    - add a CytobandTable loader + lookup
    - add required --cytoband-bed argument
    - off-panel breakpoint side -> cytoband label (e.g. 8q24.21), not OFF_PANEL
    - add gene_a_source / gene_b_source provenance columns
      (values: panel | cytoband | coordinate)

  bin/merge_translocations.py
    - carry the source token through canonical_ends() _ends tuples
    - replace best_gene() (string test against OFF_PANEL, which no longer
      exists) with best_end() ranked by provenance: panel > cytoband >
      coordinate
    - emit gene_a_source / gene_b_source in the merged output

Behaviour:
  - Writes a timestamped .bak next to each file before editing.
  - Refuses to apply twice (detects a sentinel marker in each file).
  - Applies by exact string match; if any anchor string is not found, it
    aborts that file with a clear message and leaves it untouched.

Usage (from repo root):
  python3 apply_cytoband_partner_annotation.py \\
      --repo /goast/mm-awgs-nextflow

Then, separately (this script does NOT edit configs):
  - add params.cytoband_bed_t2t to nextflow.config
  - add --cytoband-bed to the annotate module .nf
  (both printed as reminders at the end)
"""
import argparse
import datetime as _dt
import os
import shutil
import sys

SENTINEL = "# [cytoband-partner-annotation applied]"


def _read(path):
    with open(path, "r") as fh:
        return fh.read()


def _write(path, text):
    with open(path, "w") as fh:
        fh.write(text)


def _backup(path):
    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.bak_cytoband_{stamp}"
    shutil.copy2(path, bak)
    return bak


def _require(text, anchor, path):
    if anchor not in text:
        sys.stderr.write(
            f"ERROR: expected anchor not found in {path}:\n  {anchor!r}\n"
            f"File may already differ from the reviewed version. Aborting "
            f"this file untouched.\n")
        return False
    return True


# ---------------------------------------------------------------------------
# Blocks inserted into annotate_mm_translocations.py
# ---------------------------------------------------------------------------

ANNOTATE_CYTOBAND_CLASS = '''
@dataclass
class CytobandTable:
    """Cytoband lookup for T2T-CHM13v2.0. Bands tile each chromosome
    contiguously, so any in-range coordinate resolves to exactly one band."""
    # {chrom: [(start, end, band), ...] sorted by start}
    bands: Dict[str, List[Tuple[int, int, str]]]

    def band_for(self, chrom: Optional[str], pos: Optional[int]) -> Optional[str]:
        if chrom is None or pos is None:
            return None
        arr = self.bands.get(chrom)
        if not arr:
            return None
        # Linear scan is fine: <= ~60 bands per chromosome.
        for start, end, name in arr:
            if start <= pos < end:
                return name
        return None


def load_cytobands(bed_path: Path) -> CytobandTable:
    """Load a 5-column UCSC cytoband BED (chrom start end band gieStain).
    Only the first four columns are used."""
    bands: Dict[str, List[Tuple[int, int, str]]] = {}
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\\t")
            if len(parts) < 4:
                continue
            chrom, start, end, band = parts[0], int(parts[1]), int(parts[2]), parts[3]
            bands.setdefault(chrom, []).append((start, end, band))
    for chrom in bands:
        bands[chrom].sort(key=lambda t: t[0])
    if not bands:
        sys.stderr.write(f"ERROR: no cytobands parsed from {bed_path}\\n")
        sys.exit(1)
    return CytobandTable(bands)


def _strip_chr(chrom: str) -> str:
    """chr8 -> 8, chrX -> X, for clinical band notation (8q24.21)."""
    return chrom[3:] if chrom.startswith("chr") else chrom


def characterize_side(chrom, pos, region, cytobands):
    """Return (label, source) for one breakpoint side.

    panel region      -> (region.name, "panel")
    else band found   -> ("8q24.21",   "cytoband")
    else (defensive)  -> ("chr8:127.0Mb", "coordinate")
    """
    if region is not None:
        return region.name, "panel"
    band = cytobands.band_for(chrom, pos)
    if band is not None:
        return f"{_strip_chr(chrom)}{band}", "cytoband"
    if chrom is not None and pos is not None:
        return f"{chrom}:{pos / 1e6:.1f}Mb", "coordinate"
    return "OFF_PANEL", "coordinate"
'''

# Two smaller anchors (the real file has blank lines between sub-blocks, so a
# single large block match is fragile). Match the def signature and the two
# gene assignment lines independently.
ANNOTATE_SIG_OLD = '''def annotate(records, panel, dictionary, sample):'''
ANNOTATE_SIG_NEW = '''def annotate(records, panel, dictionary, sample, cytobands):'''

ANNOTATE_ANNOTATE_BODY_OLD = '''        gene_a = side_a.name if side_a else "OFF_PANEL"
        gene_b = side_b.name if side_b else "OFF_PANEL"'''

ANNOTATE_ANNOTATE_BODY_NEW = '''        gene_a, gene_a_source = characterize_side(r.chrom, r.pos, side_a, cytobands)
        gene_b, gene_b_source = characterize_side(r.mate_chrom, r.mate_pos, side_b, cytobands)'''

# The dictionary key must only match when BOTH sides are curated panel genes;
# a cytoband label must never key into the known-pair dictionary.
ANNOTATE_DICT_GUARD_OLD = '''        if side_a and side_b:
            key = tuple(sorted([gene_a.upper(), gene_b.upper()]))'''

ANNOTATE_DICT_GUARD_NEW = '''        if side_a and side_b and gene_a_source == "panel" and gene_b_source == "panel":
            key = tuple(sorted([gene_a.upper(), gene_b.upper()]))'''

ANNOTATE_ROW_OLD = '''            "gene_b":         gene_b,
            "known_mm_pair":  known,'''

ANNOTATE_ROW_NEW = '''            "gene_b":         gene_b,
            "gene_a_source":  gene_a_source,
            "gene_b_source":  gene_b_source,
            "known_mm_pair":  known,'''

ANNOTATE_CALL_OLD = '''    records = parse_vcf(args.vcf)
    rows = annotate(records, panel, dictionary, args.sample)'''

ANNOTATE_CALL_NEW = '''    records = parse_vcf(args.vcf)
    cytobands = load_cytobands(args.cytoband_bed)
    rows = annotate(records, panel, dictionary, args.sample, cytobands)'''

ANNOTATE_ARG_OLD = '''    ap.add_argument("--dictionary",  required=True, type=Path)'''

ANNOTATE_ARG_NEW = '''    ap.add_argument("--dictionary",  required=True, type=Path)
    ap.add_argument("--cytoband-bed", required=True, type=Path,
                    help="T2T-CHM13v2.0 cytoband BED (chrom start end band ...). "
                         "Off-panel breakpoint partners are characterized by band.")'''

ANNOTATE_COLS_OLD = '''        "chrom_b", "pos_b", "gene_b",
        "known_mm_pair", "known_freq",'''

ANNOTATE_COLS_NEW = '''        "chrom_b", "pos_b", "gene_b",
        "gene_a_source", "gene_b_source",
        "known_mm_pair", "known_freq",'''


def patch_annotate(repo):
    path = os.path.join(repo, "bin", "annotate_mm_translocations.py")
    if not os.path.isfile(path):
        sys.stderr.write(f"ERROR: not found: {path}\n")
        return False
    text = _read(path)
    if SENTINEL in text:
        sys.stderr.write(f"SKIP (already applied): {path}\n")
        return True

    anchors = [
        ANNOTATE_SIG_OLD, ANNOTATE_ANNOTATE_BODY_OLD, ANNOTATE_DICT_GUARD_OLD,
        ANNOTATE_ROW_OLD, ANNOTATE_CALL_OLD, ANNOTATE_ARG_OLD,
        ANNOTATE_COLS_OLD,
        "def region_for(", "def main() -> int:",
    ]
    for a in anchors:
        if not _require(text, a, path):
            return False

    bak = _backup(path)

    # Insert cytoband class/loader just before region_for().
    text = text.replace("def region_for(",
                        ANNOTATE_CYTOBAND_CLASS + "\n\ndef region_for(", 1)
    text = text.replace(ANNOTATE_SIG_OLD, ANNOTATE_SIG_NEW, 1)
    text = text.replace(ANNOTATE_ANNOTATE_BODY_OLD, ANNOTATE_ANNOTATE_BODY_NEW, 1)
    text = text.replace(ANNOTATE_DICT_GUARD_OLD, ANNOTATE_DICT_GUARD_NEW, 1)
    text = text.replace(ANNOTATE_ROW_OLD, ANNOTATE_ROW_NEW, 1)
    text = text.replace(ANNOTATE_CALL_OLD, ANNOTATE_CALL_NEW, 1)
    text = text.replace(ANNOTATE_ARG_OLD, ANNOTATE_ARG_NEW, 1)
    text = text.replace(ANNOTATE_COLS_OLD, ANNOTATE_COLS_NEW, 1)
    text = text.replace('__version__ = "0.1.0"',
                        '__version__ = "0.2.0"\n' + SENTINEL, 1)

    _write(path, text)
    sys.stderr.write(f"PATCHED: {path}\n  backup: {bak}\n")
    return True


# ---------------------------------------------------------------------------
# merge_translocations.py
# ---------------------------------------------------------------------------

MERGE_ENDS_OLD = '''    a = (row["chrom_a"], to_int(row["pos_a"], 0), row.get("gene_a", ""))
    b = (row["chrom_b"], to_int(row["pos_b"], 0), row.get("gene_b", ""))'''

MERGE_ENDS_NEW = '''    a = (row["chrom_a"], to_int(row["pos_a"], 0), row.get("gene_a", ""),
         (row.get("gene_a_source", "") or "coordinate"))
    b = (row["chrom_b"], to_int(row["pos_b"], 0), row.get("gene_b", ""),
         (row.get("gene_b_source", "") or "coordinate"))'''

MERGE_BESTGENE_OLD = '''def best_gene(members, end_index):
    """Pick a gene for end 0 or 1: prefer a non-empty, non-OFF_PANEL value."""
    fallback = ""
    for m in members:
        g = m["_ends"][end_index][2].strip()
        if g and g != "OFF_PANEL":
            return g
        if g and not fallback:
            fallback = g
    return fallback'''

MERGE_BESTGENE_NEW = '''_SOURCE_RANK = {"panel": 0, "cytoband": 1, "coordinate": 2}


def best_end(members, end_index):
    """Pick the best (gene, source) for end 0 or 1 across clustered members.

    Prefers the strongest provenance (panel > cytoband > coordinate). This
    replaces the old string test against "OFF_PANEL": off-panel partners are
    now labelled by cytoband, so the merge must rank on the provenance token
    rather than a magic gene name. Ties are broken by taking the end from the
    representative member (most callers / most support), so an adjacent-band
    disagreement resolves to the highest-support call rather than input order.
    """
    best = None
    best_rank = None
    ranked = sorted(members, key=_rep_key, reverse=True)
    for m in ranked:
        end = m["_ends"][end_index]
        gene = end[2].strip()
        source = (end[3] if len(end) > 3 else "coordinate").strip() or "coordinate"
        if not gene:
            continue
        rank = _SOURCE_RANK.get(source, 3)
        if best is None or rank < best_rank:
            best, best_rank = (gene, source), rank
    return best if best is not None else ("", "coordinate")


def best_gene(members, end_index):
    """Backwards-compatible shim: return only the gene label."""
    return best_end(members, end_index)[0]'''

MERGE_CLUSTEROUT_OLD = '''        "chrom_a": end1[0], "pos_a": str(end1[1]), "gene_a": best_gene(members, 0),
        "chrom_b": end2[0], "pos_b": str(end2[1]), "gene_b": best_gene(members, 1),'''

# Call best_end inline. It is cheap (a short loop over cluster members), so
# calling it twice per end to fetch gene then source avoids needing a separate
# _ga/_gb prep line and its fragile block anchor.
MERGE_CLUSTEROUT_NEW = '''        "chrom_a": end1[0], "pos_a": str(end1[1]),
        "gene_a": best_end(members, 0)[0],
        "gene_a_source": best_end(members, 0)[1],
        "chrom_b": end2[0], "pos_b": str(end2[1]),
        "gene_b": best_end(members, 1)[0],
        "gene_b_source": best_end(members, 1)[1],'''


def patch_merge(repo):
    path = os.path.join(repo, "bin", "merge_translocations.py")
    if not os.path.isfile(path):
        sys.stderr.write(f"ERROR: not found: {path}\n")
        return False
    text = _read(path)
    if SENTINEL in text:
        sys.stderr.write(f"SKIP (already applied): {path}\n")
        return True

    anchors = [MERGE_ENDS_OLD, MERGE_BESTGENE_OLD, MERGE_CLUSTEROUT_OLD]
    for a in anchors:
        if not _require(text, a, path):
            return False

    bak = _backup(path)
    text = text.replace(MERGE_ENDS_OLD, MERGE_ENDS_NEW, 1)
    text = text.replace(MERGE_BESTGENE_OLD, MERGE_BESTGENE_NEW, 1)
    text = text.replace(MERGE_CLUSTEROUT_OLD, MERGE_CLUSTEROUT_NEW, 1)
    # Mark applied (find a safe top-of-file spot after the module docstring).
    text = text.replace("\nimport ", "\n" + SENTINEL + "\nimport ", 1)

    _write(path, text)
    sys.stderr.write(f"PATCHED: {path}\n  backup: {bak}\n")
    sys.stderr.write(
        "  NOTE: merged output now carries gene_a_source / gene_b_source. "
        "Ensure out_cols in this file includes them (see reminder below).\n")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True,
                    help="Path to mm-awgs-nextflow repo root.")
    args = ap.parse_args()

    ok_a = patch_annotate(args.repo)
    ok_m = patch_merge(args.repo)

    sys.stderr.write("\n" + "=" * 70 + "\n")
    sys.stderr.write("MANUAL STEPS NOT DONE BY THIS SCRIPT:\n\n")
    sys.stderr.write(
        "1. merge_translocations.py out_cols: NO ACTION NEEDED. out_cols is\n"
        "   built as in_cols + [...], and in_cols is the annotate output\n"
        "   header, which now carries gene_a_source / gene_b_source. The\n"
        "   provenance columns therefore reach the merged output and the\n"
        "   dashboard automatically.\n\n")
    sys.stderr.write(
        "2. nextflow.config params block, add:\n"
        "   cytoband_bed_t2t = \"${projectDir}/assets/"
        "chm13v2.0_cytobands_allchrs.bed\"\n\n")
    sys.stderr.write(
        "3. modules/local/annotate_mm_translocations.nf script block: add\n"
        "   --cytoband-bed ${params.cytoband_bed_t2t} \\\n"
        "   to the annotate_mm_translocations.py invocation.\n")
    sys.stderr.write("=" * 70 + "\n")

    if not (ok_a and ok_m):
        sys.exit(1)


if __name__ == "__main__":
    main()
