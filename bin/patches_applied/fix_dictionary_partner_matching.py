#!/usr/bin/env python3
"""
fix_dictionary_partner_matching.py
==================================
Fix the known_mm_pair / known_freq "always empty" bug in
bin/annotate_mm_translocations.py.

Root cause: the translocation dictionary names Ig partners with a "_locus"
suffix (IGH_locus) and joins compound partners with "_" (FGFR3_NSD2,
WWOX_MAF), while the v7 panel BED labels the same regions plainly (IGH) and
with "/" (FGFR3/NSD2, WWOX/MAF). The lookup compared these as exact strings,
so IGH != IGH_locus and every canonical Ig pair missed -> known_mm_pair blank
across all calls.

Fix: match partners on gene identity rather than exact label. Each partner
name is reduced to a set of normalized gene tokens (split on / and _,
uppercase, drop the region-suffix word "LOCUS"). Two partners match when
their token sets INTERSECT, which also handles the case where the BED carries
an extra component the dictionary does not name (IGL/IGLL5 vs IGL_locus).

Changes, all in bin/annotate_mm_translocations.py:
  1. add _norm_tokens() helper
  2. load_dictionary(): return a list of (tokset_a, tokset_b, row) instead of
     a {sorted-tuple: row} dict
  3. annotate(): replace the exact-key lookup with an intersection scan; the
     existing both-ends-must-be-panel guard is preserved

Idempotent: writes a .bak, refuses to apply twice.

Usage (from repo root):
  python3 fix_dictionary_partner_matching.py --repo /goast/mm-awgs-nextflow
"""
import argparse
import datetime as _dt
import os
import shutil
import sys

SENTINEL = "# [dictionary-token-matching applied]"

# ---------------------------------------------------------------------------
# 1. Helper inserted just before load_dictionary.
# ---------------------------------------------------------------------------
HELPER_BLOCK = '''# Region-suffix words that are not gene symbols (e.g. the "locus" in
# "IGH_locus"). Dropped during tokenization so IGH_locus and IGH match.
_PARTNER_SUFFIX_TOKENS = {"LOCUS"}


def _norm_tokens(name: str) -> frozenset:
    """Reduce a partner label to a set of normalized gene tokens.

    Splits on '/' and '_' (and '+'), uppercases, and drops region-suffix
    words. E.g. 'IGH_locus' -> {IGH}; 'FGFR3/NSD2' -> {FGFR3, NSD2};
    'IGL/IGLL5' -> {IGL, IGLL5}. Two partners are considered the same locus
    when their token sets intersect.
    """
    raw = str(name).strip().upper().replace("+", "/").replace("_", "/")
    return frozenset(t for t in raw.split("/")
                     if t and t not in _PARTNER_SUFFIX_TOKENS)


'''

# ---------------------------------------------------------------------------
# 2. load_dictionary: dict -> list of (tokset_a, tokset_b, row)
# ---------------------------------------------------------------------------
LOAD_OLD = '''def load_dictionary(dict_path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    """
    Load the MM translocation dictionary, keyed by unordered (a, b) pair
    of uppercase symbols. Missing dictionary file is non-fatal.
    """
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    if not dict_path.exists():
        return out
    with open(dict_path) as fh:
        header = fh.readline().rstrip("\\n").split("\\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\\n").split("\\t")))
            a = (row.get("partner_a") or "").strip().upper()
            b = (row.get("partner_b") or "").strip().upper()
            if not a or not b:
                continue
            key = tuple(sorted([a, b]))
            out[key] = row
    return out'''

LOAD_NEW = '''def load_dictionary(dict_path):
    """
    Load the MM translocation dictionary as a list of
    (tokset_a, tokset_b, row) entries. Partner names are reduced to gene
    token sets (see _norm_tokens) so lookup matches on gene identity rather
    than exact label. Missing dictionary file is non-fatal (empty list).
    """
    out = []
    if not dict_path.exists():
        return out
    with open(dict_path) as fh:
        header = fh.readline().rstrip("\\n").split("\\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\\n").split("\\t")))
            a = (row.get("partner_a") or "").strip()
            b = (row.get("partner_b") or "").strip()
            if not a or not b:
                continue
            out.append((_norm_tokens(a), _norm_tokens(b), row))
    return out


def dictionary_lookup(dictionary, gene_a, gene_b):
    """Return the dictionary row for an unordered gene-a / gene-b pair, or
    None. Matches when the pair's token sets intersect the entry's token
    sets on both sides (in either orientation)."""
    ta, tb = _norm_tokens(gene_a), _norm_tokens(gene_b)
    if not ta or not tb:
        return None
    for da, db, row in dictionary:
        if (ta & da and tb & db) or (ta & db and tb & da):
            return row
    return None'''

# ---------------------------------------------------------------------------
# 3. annotate(): replace the exact-key lookup with the intersection scan.
#    The both-ends-panel guard added by the cytoband patch is preserved.
# ---------------------------------------------------------------------------
LOOKUP_OLD = '''        if side_a and side_b and gene_a_source == "panel" and gene_b_source == "panel":
            key = tuple(sorted([gene_a.upper(), gene_b.upper()]))
            hit = dictionary.get(key)
            if hit:'''

LOOKUP_NEW = '''        if side_a and side_b and gene_a_source == "panel" and gene_b_source == "panel":
            hit = dictionary_lookup(dictionary, gene_a, gene_b)
            if hit:'''


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", required=True, help="Path to mm-awgs-nextflow repo root.")
    args = ap.parse_args()

    path = os.path.join(args.repo, "bin", "annotate_mm_translocations.py")
    if not os.path.isfile(path):
        sys.stderr.write(f"ERROR: not found: {path}\n")
        sys.exit(1)

    with open(path) as fh:
        text = fh.read()

    if SENTINEL in text:
        sys.stderr.write(f"SKIP (already applied): {path}\n")
        return

    for anchor, label in [(LOAD_OLD, "load_dictionary"),
                          (LOOKUP_OLD, "annotate lookup"),
                          ("def load_dictionary", "load_dictionary def")]:
        if anchor not in text:
            sys.stderr.write(
                f"ERROR: expected anchor for {label} not found in {path}.\n"
                f"File differs from the reviewed version; not modifying.\n")
            sys.exit(1)

    stamp = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.bak_dictfix_{stamp}"
    shutil.copy2(path, bak)

    # Insert helper before load_dictionary, then swap the two blocks.
    text = text.replace("def load_dictionary",
                        HELPER_BLOCK + "def load_dictionary", 1)
    text = text.replace(LOAD_OLD, LOAD_NEW, 1)
    text = text.replace(LOOKUP_OLD, LOOKUP_NEW, 1)
    # Bump version + sentinel (version line already bumped by cytoband patch).
    text = text.replace(SENTINEL, "", 0)  # no-op safety
    text = text.replace('__version__ = "0.2.0"',
                        '__version__ = "0.3.0"\n' + SENTINEL, 1)
    if SENTINEL not in text:
        # cytoband patch may have used a different version string; append near top
        text = text.replace("import sys\n", "import sys\n" + SENTINEL + "\n", 1)

    with open(path, "w") as fh:
        fh.write(text)

    sys.stderr.write(f"PATCHED: {path}\n  backup: {bak}\n")
    sys.stderr.write("  known_mm_pair lookup now matches on gene tokens "
                     "(IGH_locus <-> IGH, FGFR3_NSD2 <-> FGFR3/NSD2, etc).\n")


if __name__ == "__main__":
    main()
