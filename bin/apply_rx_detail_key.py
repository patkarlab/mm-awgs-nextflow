#!/usr/bin/env python3
"""
apply_rx_detail_key.py  (rx_detail_key_v1)

Explain the fields in a rearrangement locus's detail row.

Clicking a locus opens a two-column table of raw fields from the
translocations TSV (sv_id, merged_sv_ids, supp_vec, match_quality,
anchor_class, gene_a_dist, ...) with no explanation. Each field name now
carries its explanation as a hover title with a dotted underline, a
collapsible "Detail field key" joins the Column key above the table, and the
footnote under the details is rewritten in plain terms.

Files: parsers/rearrangements.py (DETAIL_KEY, passed as `detail_key`),
templates/rearrangements_tab.html.j2. Backups <file>.dk.bak.

Usage (from the repo root):
  python3 bin/apply_rx_detail_key.py            apply
  python3 bin/apply_rx_detail_key.py --check
  python3 bin/apply_rx_detail_key.py --revert
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "rx_detail_key_v1"
BAK_SUFFIX = ".dk.bak"
DASH = "bin/dashboard_builder/"

DETAIL_KEY_SRC = '''# rx_detail_key_v1: what each field in a locus's detail row means. Order is
# the display order; fields absent from a row are simply not shown.
DETAIL_KEY = [
    ("sv_id", "identifier of the representative record, the caller-named breakend with the most supporting reads among those merged into this junction"),
    ("merged_sv_ids", "every caller record SURVIVOR merged into this junction, by the caller's own name for it; one caller can contribute two records when it wrote both mates"),
    ("n_merged", "how many caller records were merged"),
    ("support_representative", "supporting reads of the representative record, as its caller reported them"),
    ("support_sniffles", "supporting reads Sniffles reported for this junction (INFO/SUPPORT); blank when Sniffles did not call it"),
    ("support_cutesv", "supporting reads CuteSV reported (INFO/RE)"),
    ("support_severus", "supporting reads Severus reported (FORMAT/DV)"),
    ("support_savana", "supporting reads SAVANA reported (INFO/TUMOUR_READ_SUPPORT)"),
    ("filter_sniffles", "Sniffles' FILTER for its record: PASS, or GT when it could not assign a genotype, which every low-support call receives"),
    ("filter_cutesv", "CuteSV's FILTER: PASS, or q5 when its quality score is below 5, which every two- or three-read call receives"),
    ("filter_severus", "Severus' FILTER for its record"),
    ("filter_savana", "SAVANA's FILTER: PASS, or PREDICTED_NOISE when its classifier scored the call as noise"),
    ("support_nanomonsv", "supporting reads nanomonsv found for the same junction, in its independent run against the HPRC control panel"),
    ("filter_nanomonsv", "nanomonsv's FILTER: PASS means confirmed; anything else names why it rejected the junction"),
    ("supp_vec", "SURVIVOR's presence vector at merge time, one digit per caller in the order they were merged (Sniffles, CuteSV, Severus, SAVANA). The per-caller support fields above are re-derived afterwards by matching each caller's own VCF at the position, so a caller can show reads there while its digit here is 0; trust the support fields"),
    ("pon_freq", "how often a junction at this position appears in the panel of normal genomes; a common one is a mapping artefact or a germline variant, not a tumour event"),
    ("known_mm_pair", "the named myeloma rearrangement this junction matches in the reference dictionary, in ISCN form"),
    ("entity", "the WHO 5th edition entity that the named pair defines, where it defines one"),
    ("tier", "dictionary tier of the named pair: defining (a primary, disease-defining translocation), secondary (a progression event such as a MYC rearrangement), or lower"),
    ("known_freq", "published frequency of the named pair in newly diagnosed myeloma"),
    ("match_quality", "how the junction matched the dictionary entry: full when both breakends fall in the expected windows; partial when only one does and the other is nearby"),
    ("anchor", "the locus this junction is grouped under: an immunoglobulin locus (IGH_locus, IGK_locus, IGL_locus) or MYC"),
    ("anchor_class", "mechanism class of the anchor: enhancer_hijack for an immunoglobulin locus placing its enhancer next to the partner gene; for MYC the class describes the partner (Ig or non-Ig) rather than a mechanism"),
    ("flag", "why the junction is not reportable; hover the flag badge in the table for the same text"),
    ("dict_notes", "free-text note from the reference dictionary entry"),
    ("ig_region_a", "for a breakend in an immunoglobulin locus, which sub-region it falls in: C_switch (constant/switch, where primary translocations arise), J, D or V"),
    ("ig_region_b", "as ig_region_a, for the other breakend"),
    ("gene_a_dist", "distance in bp from breakend A to the gene it is named after; 0 means inside the gene body"),
    ("gene_b_dist", "as gene_a_dist, for breakend B"),
    ("filter", "the merged record's own FILTER after the ensemble, PASS unless every caller rejected it"),
]
DETAIL_KEY_MAP = dict(DETAIL_KEY)

'''

EDITS = [
    (
        DASH + "parsers/rearrangements.py",
        [
            ("COLUMN_KEY = [\n", DETAIL_KEY_SRC + "COLUMN_KEY = [\n"),
            (
                "                \"reading_key\": READING_KEY, \"column_key\": COLUMN_KEY,\n"
                "                \"other_svs\": [], \"locus_tol\": locus_tol}\n",
                "                \"reading_key\": READING_KEY, \"column_key\": COLUMN_KEY,\n"
                "                \"detail_key\": DETAIL_KEY, \"detail_key_map\": DETAIL_KEY_MAP,\n"
                "                \"other_svs\": [], \"locus_tol\": locus_tol}\n",
            ),
            (
                "        \"reading_key\": READING_KEY, \"column_key\": COLUMN_KEY,\n"
                "        \"locus_tol\": locus_tol,\n",
                "        \"reading_key\": READING_KEY, \"column_key\": COLUMN_KEY,\n"
                "        \"detail_key\": DETAIL_KEY, \"detail_key_map\": DETAIL_KEY_MAP,\n"
                "        \"locus_tol\": locus_tol,\n",
            ),
        ],
    ),
    (
        DASH + "templates/rearrangements_tab.html.j2",
        [
            (
                "        {% for head, text in rx.column_key %}\n"
                "          <tr><td class=\"text-nowrap\"><code>{{ head }}</code></td><td>{{ text }}</td></tr>\n"
                "        {% endfor %}\n"
                "      </tbody>\n"
                "    </table>\n"
                "  </details>\n",
                "        {% for head, text in rx.column_key %}\n"
                "          <tr><td class=\"text-nowrap\"><code>{{ head }}</code></td><td>{{ text }}</td></tr>\n"
                "        {% endfor %}\n"
                "      </tbody>\n"
                "    </table>\n"
                "  </details>\n"
                "  {# rx_detail_key_v1 #}\n"
                "  <details class=\"mb-3\">\n"
                "    <summary class=\"small\">Detail field key &mdash; the fields shown when a locus row is opened</summary>\n"
                "    <table class=\"table table-sm table-borderless small mb-0\" style=\"max-width:1100px\">\n"
                "      <tbody>\n"
                "        {% for head, text in rx.detail_key %}\n"
                "          <tr><td class=\"text-nowrap\"><code>{{ head }}</code></td><td>{{ text }}</td></tr>\n"
                "        {% endfor %}\n"
                "      </tbody>\n"
                "    </table>\n"
                "  </details>\n",
            ),
            (
                "              {% for k, v in l.detail %}\n"
                "                <tr><td class=\"text-muted text-nowrap pe-3\">{{ k }}</td><td style=\"overflow-wrap:anywhere\">{{ v }}</td></tr>\n",
                "              {% for k, v in l.detail %}\n"
                "                <tr><td class=\"text-muted text-nowrap pe-3\"{% if rx.detail_key_map and k in rx.detail_key_map %} title=\"{{ rx.detail_key_map[k] }}\" style=\"cursor:help;text-decoration:underline dotted\"{% endif %}>{{ k }}</td><td style=\"overflow-wrap:anywhere\">{{ v }}</td></tr>\n",
            ),
            (
                "            <p class=\"mb-0\">Fields are those of the best-supported breakend; per-caller read counts are re-derived from each caller's own VCF.</p>\n",
                "            <p class=\"mb-0\">These fields belong to the breakend with the most supporting reads in this locus. Read counts per caller are taken from each caller's own output, not from the merged record. Hover a field name for its meaning, or open the Detail field key above the table.</p>\n",
            ),
        ],
    ),
]


def state(root):
    out = []
    for rel, edits in EDITS:
        p = root / rel
        if not p.is_file():
            out.append((rel, "MISSING FILE")); continue
        text = p.read_text()
        if SENTINEL in text:
            out.append((rel, "applied")); continue
        bad = [a for a, _ in edits if text.count(a) != 1]
        if bad:
            out.append((rel, "ANCHOR MISMATCH: %d of %d anchors not unique" % (len(bad), len(edits))))
            for a in bad:
                sys.stderr.write("    count=%d  %r\n" % (text.count(a), a.splitlines()[0]))
            continue
        out.append((rel, "ready"))
    return out


def apply(root):
    st = state(root)
    if any(s not in ("ready", "applied") for _, s in st):
        for rel, s in st:
            print("  %-52s %s" % (rel, s))
        sys.exit("Refusing to apply: fix the problems above first.")
    if all(s == "applied" for _, s in st):
        print("Already applied; nothing to do."); return
    for rel, edits in EDITS:
        p = root / rel
        text = p.read_text()
        if SENTINEL in text:
            print("  %-52s already applied, skipped" % rel); continue
        bak = p.with_name(p.name + BAK_SUFFIX)
        shutil.copy2(p, bak)
        for a, r in edits:
            text = text.replace(a, r, 1)
        assert SENTINEL in text
        p.write_text(text)
        print("  %-52s patched (backup: %s)" % (rel, bak.name))


def revert(root):
    for rel, _ in EDITS:
        p = root / rel
        bak = p.with_name(p.name + BAK_SUFFIX)
        if not bak.is_file():
            print("  %-52s no %s, skipped" % (rel, BAK_SUFFIX)); continue
        shutil.copy2(bak, p); bak.unlink()
        print("  %-52s restored from %s" % (rel, BAK_SUFFIX))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1].strip())
    ap.add_argument("--root", default=".")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true")
    g.add_argument("--revert", action="store_true")
    args = ap.parse_args()
    root = Path(args.root).resolve()
    if args.check:
        for rel, s in state(root):
            print("  %-52s %s" % (rel, s))
        return
    if args.revert:
        revert(root); return
    apply(root)


if __name__ == "__main__":
    main()
