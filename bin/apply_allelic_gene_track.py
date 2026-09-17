#!/usr/bin/env python3
"""
apply_allelic_gene_track.py  (allelic_gene_track_v1)

Gene boundaries and a legend on the per-window depth + allele-fraction plot
of the "Copy number & allelic state" tab.

What changes
------------
The plot for an assessed window (TP53+TNFSF12, WWOX/MAF, MYC, ...) gains a
gene track between the depth panel and the allele-fraction panel, drawn from
the hg38 gene model (assets/aWGS_PCN_v7_gene_model_hg38.bed): one labelled
band per panel gene overlapping the window, so a reader can see where TP53
sits inside the 995 kb window rather than infer it. The prose caption under
the plot becomes a swatch legend plus a four-line reading guide.

How the gene model reaches the page
-----------------------------------
The report is self-contained, so the model travels inside the bundle:

  workflows/mm_awgs.nf         passes params.gene_model_hg38 into REPORT_BUNDLE
                               (an empty list when the param is null)
  modules/local/report_bundle.nf  stages it and exports GENE_MODEL_HG38, the
                               same pattern as PANEL_BED_T2T / PANEL_BED_HG38
  bin/build_report_bundle.sh   copies it to <bundle>/<sample>/allelic/
                               gene_model_hg38.bed beside the genebaf files
                               (3 KB per sample); run by hand it falls back to
                               ../assets/ beside the script
  parsers/allelic.py           reads that file and passes the genes that
                               overlap each window to _window_svg
  copy_number_tab.html.j2      legend block replaces the caption

Usage (from the repo root):
  python3 bin/apply_allelic_gene_track.py            apply
  python3 bin/apply_allelic_gene_track.py --check    report state, change nothing
  python3 bin/apply_allelic_gene_track.py --revert   restore the backups

Backups are written as <file>.gt.bak. The sentinel `allelic_gene_track_v1`
refuses a second application. Requires bundle_cn_kind_v1 to be present in
build_report_bundle.sh (this patch anchors below it).
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "allelic_gene_track_v1"
BAK_SUFFIX = ".gt.bak"

DASH = "bin/dashboard_builder/"

EDITS = [
    (
        "workflows/mm_awgs.nf",
        [
            (
                "            file(params.panel_bed_hg38, checkIfExists: true),\n"
                "            params.outdir,\n"
                "            bundle_name\n"
                "        )\n",
                "            file(params.panel_bed_hg38, checkIfExists: true),\n"
                "            // allelic_gene_track_v1: hg38 gene model for the window plots.\n"
                "            // Optional; an empty list stages nothing and the bundle script\n"
                "            // draws no gene track.\n"
                "            params.gene_model_hg38 ? file(params.gene_model_hg38, checkIfExists: true) : [],\n"
                "            params.outdir,\n"
                "            bundle_name\n"
                "        )\n",
            ),
        ],
    ),
    (
        "modules/local/report_bundle.nf",
        [
            (
                "    path panel_bed_t2t\n"
                "    path panel_bed_hg38\n"
                "    val  outdir\n",
                "    path panel_bed_t2t\n"
                "    path panel_bed_hg38\n"
                "    // allelic_gene_track_v1: hg38 gene model, copied into each sample's\n"
                "    // allelic/ so the window plots can draw gene boundaries. May be an\n"
                "    // empty list when params.gene_model_hg38 is null.\n"
                "    path gene_model_hg38\n"
                "    val  outdir\n",
            ),
            (
                "    export PANEL_BED_HG38=\"\\$(readlink -f ${panel_bed_hg38})\"\n",
                "    export PANEL_BED_HG38=\"\\$(readlink -f ${panel_bed_hg38})\"\n"
                "    export GENE_MODEL_HG38=\"${gene_model_hg38 ? '\\$(readlink -f ' + gene_model_hg38 + ')' : ''}\"\n",
            ),
        ],
    ),
    (
        "bin/build_report_bundle.sh",
        [
            (
                "RESULTS=\"${1:?Usage: build_report_bundle.sh <results_dir> [bundle_name]}\"\n",
                "# allelic_gene_track_v1: hg38 gene model for the allelic window plots.\n"
                "# Inside a Nextflow task it is staged and exported by REPORT_BUNDLE; run\n"
                "# by hand it sits beside this script under ../assets. Empty or missing\n"
                "# means no gene track, not a failure.\n"
                "GENE_MODEL_HG38=\"${GENE_MODEL_HG38-${SCRIPT_DIR}/../assets/aWGS_PCN_v7_gene_model_hg38.bed}\"\n"
                "\n"
                "RESULTS=\"${1:?Usage: build_report_bundle.sh <results_dir> [bundle_name]}\"\n",
            ),
            (
                "  copy_all   \"$d/allelic\" \"$s\" -path '*hg38/allelic*' -name \"${s}.genebaf.sites.tsv\"\n",
                "  copy_all   \"$d/allelic\" \"$s\" -path '*hg38/allelic*' -name \"${s}.genebaf.sites.tsv\"\n"
                "  # allelic_gene_track_v1: the gene model rides beside the genebaf files\n"
                "  # (3 KB) so parsers/allelic.py reads only its own subdirectory.\n"
                "  if [[ -d \"$d/allelic\" && -n \"$GENE_MODEL_HG38\" && -f \"$GENE_MODEL_HG38\" ]]; then\n"
                "    cp -L \"$GENE_MODEL_HG38\" \"$d/allelic/gene_model_hg38.bed\"\n"
                "    echo \"  + allelic/gene_model_hg38.bed\"\n"
                "  fi\n",
            ),
        ],
    ),
    (
        DASH + "parsers/allelic.py",
        [
            (
                "def _window_svg(gene, start, end, bins, sites, width=960):\n"
                "    \"\"\"Inline SVG: on-target depth log2 per 1 kb bin (top) and per-site\n"
                "    alternate allele fraction (bottom) across one panel window. No\n"
                "    dependencies, self-contained in the page. (allelic_plot_v1)\"\"\"\n"
                "    if not bins and not sites:\n"
                "        return \"\"\n"
                "    span = max(1, end - start)\n"
                "    left, right, h_top, h_bot, gap, pad = 52, 12, 110, 110, 26, 16\n",
                "GENE_MODEL_NAME = \"gene_model_hg38.bed\"\n"
                "\n"
                "\n"
                "def _read_gene_model(d):\n"
                "    \"\"\"Panel gene bodies in hg38 coordinates (allelic_gene_track_v1).\n"
                "    Copied into allelic/ by build_report_bundle.sh; absent in bundles built\n"
                "    before it, in which case the plots simply carry no gene track.\"\"\"\n"
                "    path = os.path.join(d, GENE_MODEL_NAME)\n"
                "    genes = []\n"
                "    if not os.path.isfile(path):\n"
                "        return genes\n"
                "    with open(path) as fh:\n"
                "        for line in fh:\n"
                "            if not line.strip() or line.startswith((\"#\", \"track\", \"browser\")):\n"
                "                continue\n"
                "            f = line.rstrip(\"\\n\").split(\"\\t\")\n"
                "            if len(f) < 4:\n"
                "                continue\n"
                "            s, e = _int(f[1]), _int(f[2])\n"
                "            if s is None or e is None:\n"
                "                continue\n"
                "            genes.append({\"chrom\": f[0], \"start\": s, \"end\": e, \"name\": f[3]})\n"
                "    return genes\n"
                "\n"
                "\n"
                "def _genes_in_window(genes, chrom, start, end):\n"
                "    return [g for g in genes\n"
                "            if g[\"chrom\"] == chrom and g[\"end\"] > start and g[\"start\"] < end]\n"
                "\n"
                "\n"
                "def _window_svg(gene, start, end, bins, sites, genes=None, width=960):\n"
                "    \"\"\"Inline SVG: on-target depth log2 per 1 kb bin (top), a gene track\n"
                "    (middle, allelic_gene_track_v1) and per-site alternate allele fraction\n"
                "    (bottom) across one panel window. No dependencies, self-contained in\n"
                "    the page. (allelic_plot_v1)\"\"\"\n"
                "    if not bins and not sites:\n"
                "        return \"\"\n"
                "    genes = genes or []\n"
                "    span = max(1, end - start)\n"
                "    # The gap between panels holds the gene track when there is one.\n"
                "    left, right, h_top, h_bot, pad = 52, 12, 110, 110, 16\n"
                "    gap = 52 if genes else 26\n",
            ),
            (
                "    out.append('<text x=\"%d\" y=\"%.1f\" fill=\"#444\" font-weight=\"600\">alternate allele fraction at common SNPs</text>' % (left, pad + h_top + gap - 6))\n",
                "    out.append('<text x=\"%d\" y=\"%.1f\" fill=\"#444\" font-weight=\"600\">alternate allele fraction at common SNPs</text>' % (left, pad + h_top + gap - 6))\n"
                "    # gene track: one band per panel gene overlapping the window, clipped\n"
                "    # to the window, labelled at the centre of the visible part.\n"
                "    if genes:\n"
                "        y_band = pad + h_top + 8\n"
                "        out.append('<text x=\"%d\" y=\"%.1f\" text-anchor=\"end\" fill=\"#666\">genes</text>' % (left - 6, y_band + 8))\n"
                "        for g in sorted(genes, key=lambda r: r[\"start\"]):\n"
                "            gs, ge = max(g[\"start\"], start), min(g[\"end\"], end)\n"
                "            x0, x1 = x(gs), x(ge)\n"
                "            out.append('<rect x=\"%.1f\" y=\"%.1f\" width=\"%.2f\" height=\"9\" fill=\"#2e7d32\" opacity=\"0.85\"><title>%s chr:%s-%s</title></rect>'\n"
                "                       % (x0, y_band, max(1.5, x1 - x0), g[\"name\"], format(g[\"start\"], \",\"), format(g[\"end\"], \",\")))\n"
                "            out.append('<text x=\"%.1f\" y=\"%.1f\" text-anchor=\"middle\" fill=\"#1b5e20\" font-weight=\"600\">%s</text>'\n"
                "                       % ((x0 + x1) / 2.0, y_band + 21, g[\"name\"]))\n",
            ),
            (
                "    bins = _read_tsv(os.path.join(d, \"%s.genebaf.bins.tsv\" % sample))\n",
                "    bins = _read_tsv(os.path.join(d, \"%s.genebaf.bins.tsv\" % sample))\n"
                "    gene_model = _read_gene_model(d)\n",
            ),
            (
                "        row[\"svg\"] = (_window_svg(gene, row[\"start\"] or 0, row[\"end\"] or 0,\n"
                "                                  bins_by_gene.get(gene, []), sites_by_gene.get(gene, []))\n"
                "                      if row[\"baf_scope\"] == \"assessed\" else \"\")\n",
                "        row[\"genes_in_window\"] = _genes_in_window(\n"
                "            gene_model, row[\"chrom\"], row[\"start\"] or 0, row[\"end\"] or 0)\n"
                "        row[\"svg\"] = (_window_svg(gene, row[\"start\"] or 0, row[\"end\"] or 0,\n"
                "                                  bins_by_gene.get(gene, []), sites_by_gene.get(gene, []),\n"
                "                                  row[\"genes_in_window\"])\n"
                "                      if row[\"baf_scope\"] == \"assessed\" else \"\")\n",
            ),
        ],
    ),
    (
        DASH + "templates/copy_number_tab.html.j2",
        [
            (
                "                        <div class=\"small text-muted\">Red points: allele fractions away from 0 and 1 (heterozygous or imbalanced sites); grey: near-homozygous. A diploid balanced window puts the red band at 0.5; a deletion or copy-neutral LOH splits it; a gain shifts it to one side.</div>\n",
                "                        {# allelic_gene_track_v1: swatch legend and reading guide #}\n"
                "                        <div class=\"small text-muted mt-1\">\n"
                "                          <span class=\"me-3\"><span style=\"display:inline-block;width:14px;height:3px;background:#1F3B73;vertical-align:middle\"></span> depth log2, 1 kb bin</span>\n"
                "                          {% if w.genes_in_window %}<span class=\"me-3\"><span style=\"display:inline-block;width:14px;height:9px;background:#2e7d32;vertical-align:middle\"></span> gene body</span>{% endif %}\n"
                "                          <span class=\"me-3\"><span style=\"display:inline-block;width:9px;height:9px;border-radius:50%;background:#b91c1c;vertical-align:middle\"></span> heterozygous SNP (allele fraction 0.15&ndash;0.85)</span>\n"
                "                          <span class=\"me-3\"><span style=\"display:inline-block;width:9px;height:9px;border-radius:50%;background:#888;vertical-align:middle\"></span> homozygous SNP (near 0 or 1; carries no allelic information)</span>\n"
                "                        </div>\n"
                "                        <div class=\"small text-muted\">\n"
                "                          Only the red sites report allelic state. <strong>Balanced diploid:</strong> one red band at 0.5, depth flat.\n"
                "                          <strong>Deletion:</strong> band splits symmetrically about 0.5, depth falls.\n"
                "                          <strong>Copy-neutral LOH:</strong> band splits, depth flat.\n"
                "                          <strong>Gain:</strong> band splits to 1/3 and 2/3, depth rises.\n"
                "                          Separation grows with the fraction of cells carrying the event; near 100% the red sites reach the grey lines.\n"
                "                          Runs of grey-only sites are haplotype blocks, not LOH.\n"
                "                        </div>\n",
            ),
        ],
    ),
]


def state(root):
    out = []
    for rel, edits in EDITS:
        p = root / rel
        if not p.is_file():
            out.append((rel, "MISSING FILE"))
            continue
        text = p.read_text()
        if SENTINEL in text:
            out.append((rel, "applied"))
            continue
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
            print("  %-48s %s" % (rel, s))
        sys.exit("Refusing to apply: fix the problems above first.")
    if all(s == "applied" for _, s in st):
        print("Already applied (sentinel present in every file); nothing to do.")
        return
    for rel, edits in EDITS:
        p = root / rel
        text = p.read_text()
        if SENTINEL in text:
            print("  %-48s already applied, skipped" % rel)
            continue
        bak = p.with_name(p.name + BAK_SUFFIX)
        shutil.copy2(p, bak)
        for anchor, repl in edits:
            text = text.replace(anchor, repl, 1)
        assert SENTINEL in text, "sentinel absent after patching %s" % rel
        p.write_text(text)
        print("  %-48s patched (backup: %s)" % (rel, bak.name))


def revert(root):
    for rel, _ in EDITS:
        p = root / rel
        bak = p.with_name(p.name + BAK_SUFFIX)
        if not bak.is_file():
            print("  %-48s no %s, skipped" % (rel, BAK_SUFFIX))
            continue
        shutil.copy2(bak, p)
        bak.unlink()
        print("  %-48s restored from %s" % (rel, BAK_SUFFIX))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1].strip())
    ap.add_argument("--root", default=".", help="repository root (default: cwd)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true", help="report state only")
    g.add_argument("--revert", action="store_true", help="restore backup copies")
    args = ap.parse_args()
    root = Path(args.root).resolve()

    if args.check:
        for rel, s in state(root):
            print("  %-48s %s" % (rel, s))
        return
    if args.revert:
        revert(root)
        return
    apply(root)


if __name__ == "__main__":
    main()
