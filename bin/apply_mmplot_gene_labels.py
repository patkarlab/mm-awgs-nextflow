#!/usr/bin/env python3
"""
apply_mmplot_gene_labels.py  (mmplot_gene_labels_v1)

Stagger the gene labels above the depth panel on single-chromosome pages so
neighbouring genes do not overprint.

Why: the hg38 gene model (de897cd) names genes individually, so pairs the
panel captures in one window are drawn as two bands 60-100 kb apart:
FGFR3/NSD2, TP53/TNFSF12, NRAS/TENT5C, IGL/IGLL5, FCRL5/FCRL4, and the
IKZF3-STAT3-MAP3K14 run on 17q. Their labels were written at the same
height and became unreadable. The targets track already staggers its labels
on two rows (target_labels_v1); the gene labels now do the same, on as many
rows as the collisions need, with the title pad grown to make room.

Method: labels are placed in order of position; each takes the lowest row
whose previous label ends before this one starts, using an estimated text
width in data units (characters x font size x 0.65 pt, converted through the
17-inch figure width). A short tick joins each raised label to its band.

Usage (from the repo root):
  python3 bin/apply_mmplot_gene_labels.py            apply
  python3 bin/apply_mmplot_gene_labels.py --check    report state, change nothing
  python3 bin/apply_mmplot_gene_labels.py --revert   restore the backup
"""

import argparse
import shutil
import sys
from pathlib import Path

SENTINEL = "mmplot_gene_labels_v1"
BAK_SUFFIX = ".gl.bak"

EDITS = [
    (
        "bin/mmplot.py",
        [
            (
                "def draw_panels(axes, chrom, length, bins, windows, segments, purity, ploidy,\n"
                "                max_cn, point_size, genes=None, targets=None, cytobands=None,\n"
                "                label_bands=False, compact=False, centromeres=None):\n"
                "    \"\"\"Fill the depth, BAF, copy number, targets and cytogram axes.\"\"\"\n",
                "def _label_rows(items, limit, figure_width_in=17.0, axes_fraction=0.86,\n"
                "                font_pt=8.0, pad=1.8):\n"
                "    \"\"\"Row index per label so neighbours do not overprint\n"
                "    (mmplot_gene_labels_v1).\n"
                "\n"
                "    items: [(centre_in_data_units, text), ...] in any order. Each label\n"
                "    takes the lowest row whose last label ends before this one starts.\n"
                "    Text width is estimated at 0.65 em per character (bold) and converted to\n"
                "    data units through the axis width, then padded.\n"
                "    \"\"\"\n"
                "    per_unit = figure_width_in * axes_fraction / max(limit, 1e-9)  # in/data\n"
                "    ends = []  # last occupied x per row\n"
                "    rows = {}\n"
                "    for centre, text in sorted(items, key=lambda t: t[0]):\n"
                "        width = len(text) * font_pt * 0.65 / 72.0 / per_unit * pad\n"
                "        start, end = centre - width / 2.0, centre + width / 2.0\n"
                "        for row, last_end in enumerate(ends):\n"
                "            if start > last_end:\n"
                "                ends[row] = end\n"
                "                rows[(centre, text)] = row\n"
                "                break\n"
                "        else:\n"
                "            ends.append(end)\n"
                "            rows[(centre, text)] = len(ends) - 1\n"
                "    return rows, len(ends)\n"
                "\n"
                "\n"
                "def draw_panels(axes, chrom, length, bins, windows, segments, purity, ploidy,\n"
                "                max_cn, point_size, genes=None, targets=None, cytobands=None,\n"
                "                label_bands=False, compact=False, centromeres=None):\n"
                "    \"\"\"Fill the depth, BAF, copy number, targets and cytogram axes.\n"
                "\n"
                "    Returns the number of gene-label rows drawn above the depth panel, so\n"
                "    the caller can pad the title to clear them.\"\"\"\n"
                "    label_row_count = 0\n",
            ),
            (
                "    # gene highlight bands behind everything\n"
                "    if genes and not compact:\n"
                "        for i, (start, end, name, _extra) in enumerate(genes.get(chrom, [])):\n"
                "            colour = HIGHLIGHT[i % len(HIGHLIGHT)]\n"
                "            for axis in (ax_depth, ax_baf, ax_cn):\n"
                "                if axis is not None:\n"
                "                    axis.axvspan(start / scale, end / scale, color=colour,\n"
                "                                 alpha=0.55, zorder=0, linewidth=0)\n"
                "            if ax_depth is not None and name:\n"
                "                ax_depth.annotate(name, ((start + end) / 2 / scale, 1.02),\n"
                "                                  xycoords=(\"data\", \"axes fraction\"),\n"
                "                                  ha=\"center\", va=\"bottom\", fontsize=8,\n"
                "                                  fontweight=\"bold\", color=\"0.2\")\n",
                "    # gene highlight bands behind everything\n"
                "    if genes and not compact:\n"
                "        gene_list = genes.get(chrom, [])\n"
                "        # mmplot_gene_labels_v1: labels staggered onto rows so genes a few\n"
                "        # tens of kb apart (TP53/TNFSF12, FGFR3/NSD2) stay legible.\n"
                "        rows, label_row_count = _label_rows(\n"
                "            [((s + e) / 2 / scale, n) for s, e, n, _x in gene_list if n], limit)\n"
                "        row_step = 0.075  # axes fraction per label row\n"
                "        for i, (start, end, name, _extra) in enumerate(gene_list):\n"
                "            colour = HIGHLIGHT[i % len(HIGHLIGHT)]\n"
                "            for axis in (ax_depth, ax_baf, ax_cn):\n"
                "                if axis is not None:\n"
                "                    axis.axvspan(start / scale, end / scale, color=colour,\n"
                "                                 alpha=0.55, zorder=0, linewidth=0)\n"
                "            if ax_depth is not None and name:\n"
                "                centre = (start + end) / 2 / scale\n"
                "                row = rows.get((centre, name), 0)\n"
                "                y = 1.02 + row * row_step\n"
                "                ax_depth.annotate(name, (centre, y),\n"
                "                                  xycoords=(\"data\", \"axes fraction\"),\n"
                "                                  ha=\"center\", va=\"bottom\", fontsize=8,\n"
                "                                  fontweight=\"bold\", color=\"0.2\",\n"
                "                                  annotation_clip=False)\n"
                "                if row > 0:\n"
                "                    # tick from the band up to a raised label\n"
                "                    ax_depth.plot([centre, centre], [1.0, y - 0.005],\n"
                "                                  transform=ax_depth.get_xaxis_transform(),\n"
                "                                  color=\"0.6\", linewidth=0.6, clip_on=False,\n"
                "                                  zorder=1)\n",
            ),
            (
                "        draw_cytogram(ax_ideo, chrom, length, cytobands or {}, label_bands,\n"
                "                      centromeres=centromeres)\n"
                "        ax_ideo.tick_params(labelbottom=True, labelsize=8)\n"
                "\n",
                "        draw_cytogram(ax_ideo, chrom, length, cytobands or {}, label_bands,\n"
                "                      centromeres=centromeres)\n"
                "        ax_ideo.tick_params(labelbottom=True, labelsize=8)\n"
                "    return label_row_count\n"
                "\n",
            ),
            (
                "        draw_panels(axes, chrom, lengths[chrom], bins, windows, segments,\n"
                "                    rho, psi, args.max_cn,\n"
                "                    point_size=13.0 * args.point_scale, genes=genes,\n"
                "                    targets=targets, cytobands=cytobands, label_bands=True,\n"
                "                    centromeres=centromeres)\n",
                "        label_rows = draw_panels(\n"
                "                    axes, chrom, lengths[chrom], bins, windows, segments,\n"
                "                    rho, psi, args.max_cn,\n"
                "                    point_size=13.0 * args.point_scale, genes=genes,\n"
                "                    targets=targets, cytobands=cytobands, label_bands=True,\n"
                "                    centromeres=centromeres) or 0\n",
            ),
            (
                "        built[0].set_title(wrap_header(\"%s  --  %s\" % (chrom, header), 17, 13),\n"
                "                           fontsize=13, pad=22)\n",
                "        # mmplot_gene_labels_v1: one extra label row needs about 11 pt.\n"
                "        built[0].set_title(wrap_header(\"%s  --  %s\" % (chrom, header), 17, 13),\n"
                "                           fontsize=13, pad=22 + 11 * max(0, label_rows - 1))\n",
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
            print("  %-20s %s" % (rel, s))
        sys.exit("Refusing to apply: fix the problems above first.")
    if all(s == "applied" for _, s in st):
        print("Already applied (sentinel present); nothing to do.")
        return
    for rel, edits in EDITS:
        p = root / rel
        text = p.read_text()
        bak = p.with_name(p.name + BAK_SUFFIX)
        shutil.copy2(p, bak)
        for anchor, repl in edits:
            text = text.replace(anchor, repl, 1)
        assert SENTINEL in text
        p.write_text(text)
        print("  %-20s patched (backup: %s)" % (rel, bak.name))


def revert(root):
    for rel, _ in EDITS:
        p = root / rel
        bak = p.with_name(p.name + BAK_SUFFIX)
        if not bak.is_file():
            print("  %-20s no %s, skipped" % (rel, BAK_SUFFIX))
            continue
        shutil.copy2(bak, p)
        bak.unlink()
        print("  %-20s restored from %s" % (rel, BAK_SUFFIX))


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
            print("  %-20s %s" % (rel, s))
        return
    if args.revert:
        revert(root)
        return
    apply(root)


if __name__ == "__main__":
    main()
