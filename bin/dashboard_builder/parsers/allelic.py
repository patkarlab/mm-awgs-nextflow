"""
allelic.py - dashboard parser joining genebaf allelic state to ichorCNA copy
number, per panel window, per chromosome.

Inputs, under <effective_dir>/allelic/:
    <sample>.genebaf.tsv        one row per panel window: depth, allelic state
    <sample>.armloh.tsv         pooled arm verdicts with covered span
    <sample>.genebaf.bins.tsv   1 kb depth bins per window (log2 vs normaliser)
    <sample>.genebaf.json       scope and parameters
plus the ichorkaryo JSON already parsed by parsers/copy_number.py, whose
`genes` list gives the ichorCNA copy number, event and cell fraction at the
same windows.

Output: {"found", "chromosomes": {chrom: {"windows": [...], "arms": [...]}},
         "scope": {...}, "reading": [...]} where each window row carries
    gene, band, size, mean_depth, log2_ontarget (median of its 1 kb bins),
    cn, cn_event, cell_fraction (from ichorCNA), baf_scope, allelic_state,
    n_het, n_usable, median_dev, expected/observed het fraction, and a
    one-word `reading` combining depth and allelic state:
        loss + LOH          -> "deletion"
        neutral + imbalanced with hets lost -> "copy-neutral LOH"
        gain + imbalanced   -> "gain (allelic)"
        neutral + balanced  -> "no event"
        anything not assessed -> depth only

Standard library only.
"""

import csv
import json
import os
import statistics
from collections import defaultdict

SUBDIR = "allelic"
CHROM_ORDER = ["chr%d" % i for i in range(1, 23)] + ["chrX", "chrY"]

READING = [
    ("Depth (on-target)", "median log2 depth ratio of the window's 1 kb bins against the sample's own diploid level, from the reads adaptive sampling kept. Independent of ichorCNA."),
    ("Copy number", "ichorCNA's integer state at the window, from off-target reads, with the cell fraction ichorkaryo derives from flow purity."),
    ("Allelic state", "at the eight wide windows only: balanced (heterozygous SNPs sit near 0.5), imbalanced (they do not), insufficient (too few usable heterozygous sites). Elsewhere: not assessed, by design — a 100 kb gene body is one or two haplotype blocks and its allele fractions cannot be read."),
    ("Reading", "depth and allelic state read together. Loss with allelic imbalance is a deletion; neutral depth with imbalance and lost heterozygosity is copy-neutral LOH; gain with imbalance is an allelic gain; neutral and balanced is no event."),
    ("Arm LOH", "the pooled heterozygosity test over every assessed window on the arm, with what those windows actually cover: n windows, the span from first to last, and the fraction of the arm inside them. One window cannot be tested against itself."),
]


def _float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _read_tsv(path):
    if not path or not os.path.isfile(path):
        return []
    with open(path, newline="") as fh:
        return [dict(r) for r in csv.DictReader(fh, delimiter="\t")]


def _reading(cn_event, state, obs_het, exp_het):
    if state not in ("balanced", "imbalanced"):
        return "depth only" if state in ("not_assessed", "", None) else "insufficient"
    lost = (obs_het is not None and exp_het and obs_het < 0.6 * exp_het)
    if cn_event == "loss":
        return "deletion" if state == "imbalanced" else "loss, alleles balanced"
    if cn_event == "gain":
        return "gain (allelic)" if state == "imbalanced" else "gain, alleles balanced"
    if state == "imbalanced":
        return "copy-neutral LOH" if lost else "allelic imbalance"
    return "no event"


GENE_MODEL_NAME = "gene_model_hg38.bed"


def _read_gene_model(d):
    """Panel gene bodies in hg38 coordinates (allelic_gene_track_v1).
    Copied into allelic/ by build_report_bundle.sh; absent in bundles built
    before it, in which case the plots simply carry no gene track."""
    path = os.path.join(d, GENE_MODEL_NAME)
    genes = []
    if not os.path.isfile(path):
        return genes
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            s, e = _int(f[1]), _int(f[2])
            if s is None or e is None:
                continue
            genes.append({"chrom": f[0], "start": s, "end": e, "name": f[3]})
    return genes


def _genes_in_window(genes, chrom, start, end):
    return [g for g in genes
            if g["chrom"] == chrom and g["end"] > start and g["start"] < end]


def _window_svg(gene, start, end, bins, sites, genes=None, width=960):
    """Inline SVG: on-target depth log2 per 1 kb bin (top), a gene track
    (middle, allelic_gene_track_v1) and per-site alternate allele fraction
    (bottom) across one panel window. No dependencies, self-contained in
    the page. (allelic_plot_v1)"""
    if not bins and not sites:
        return ""
    genes = genes or []
    span = max(1, end - start)
    # The gap between panels holds the gene track when there is one.
    left, right, h_top, h_bot, pad = 52, 12, 110, 110, 16
    gap = 52 if genes else 26
    w = width - left - right
    height = pad + h_top + gap + h_bot + 26

    def x(pos):
        return left + (pos - start) / span * w

    def y_top(v):  # log2 in [-2, 2]
        v = max(-2.0, min(2.0, v))
        return pad + (2.0 - v) / 4.0 * h_top

    def y_bot(v):  # AF in [0, 1]
        return pad + h_top + gap + (1.0 - v) * h_bot

    out = ['<svg viewBox="0 0 %d %d" width="100%%" style="max-width:%dpx;font-family:system-ui,sans-serif;font-size:10px">'
           % (width, height, width)]
    # axes and guides
    for v in (-1, 0, 1):
        out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="#ddd" stroke-dasharray="%s"/>'
                   % (left, left + w, y_top(v), y_top(v), "3,3" if v else "none"))
        out.append('<text x="%d" y="%.1f" text-anchor="end" fill="#666">%d</text>' % (left - 6, y_top(v) + 3, v))
    for v in (0.0, 0.5, 1.0):
        out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="#ddd" stroke-dasharray="%s"/>'
                   % (left, left + w, y_bot(v), y_bot(v), "3,3" if v == 0.5 else "none"))
        out.append('<text x="%d" y="%.1f" text-anchor="end" fill="#666">%.1f</text>' % (left - 6, y_bot(v) + 3, v))
    out.append('<text x="%d" y="%d" fill="#444" font-weight="600">%s &middot; depth log2 (1 kb bins)</text>' % (left, pad - 2, gene))
    out.append('<text x="%d" y="%.1f" fill="#444" font-weight="600">alternate allele fraction at common SNPs</text>' % (left, pad + h_top + gap - 6))
    # gene track: one band per panel gene overlapping the window, clipped
    # to the window, labelled at the centre of the visible part.
    if genes:
        y_band = pad + h_top + 8
        out.append('<text x="%d" y="%.1f" text-anchor="end" fill="#666">genes</text>' % (left - 6, y_band + 8))
        for g in sorted(genes, key=lambda r: r["start"]):
            gs, ge = max(g["start"], start), min(g["end"], end)
            x0, x1 = x(gs), x(ge)
            out.append('<rect x="%.1f" y="%.1f" width="%.2f" height="9" fill="#2e7d32" opacity="0.85"><title>%s chr:%s-%s</title></rect>'
                       % (x0, y_band, max(1.5, x1 - x0), g["name"], format(g["start"], ","), format(g["end"], ",")))
            out.append('<text x="%.1f" y="%.1f" text-anchor="middle" fill="#1b5e20" font-weight="600">%s</text>'
                       % ((x0 + x1) / 2.0, y_band + 21, g["name"]))
    # bins
    for b in bins:
        v = b.get("log2")
        if v is None:
            continue
        out.append('<rect x="%.1f" y="%.1f" width="%.2f" height="2" fill="#1F3B73" opacity="0.7"/>'
                   % (x(b["start"]), y_top(v) - 1, max(1.0, (b["end"] - b["start"]) / span * w)))
    # sites
    for s in sites:
        out.append('<circle cx="%.1f" cy="%.1f" r="1.8" fill="%s" opacity="0.75"/>'
                   % (x(s["pos"]), y_bot(s["af"]), "#b91c1c" if 0.15 < s["af"] < 0.85 else "#888"))
    # x axis
    for k in range(6):
        pos = start + span * k / 5
        out.append('<text x="%.1f" y="%d" text-anchor="middle" fill="#666">%.2f Mb</text>' % (x(pos), height - 8, pos / 1e6))
    out.append('</svg>')
    return "".join(out)


def parse(effective_dir, sample, ichorkaryo_genes=None):
    d = os.path.join(effective_dir, SUBDIR)
    result = {"found": False, "searched": d, "chromosomes": {}, "scope": {},
              "reading": READING, "n_assessed": 0, "n_windows": 0}
    win = _read_tsv(os.path.join(d, "%s.genebaf.tsv" % sample))
    if not win:
        return result
    arms = _read_tsv(os.path.join(d, "%s.armloh.tsv" % sample))
    bins = _read_tsv(os.path.join(d, "%s.genebaf.bins.tsv" % sample))
    gene_model = _read_gene_model(d)
    scope = {}
    jp = os.path.join(d, "%s.genebaf.json" % sample)
    if os.path.isfile(jp):
        try:
            with open(jp) as fh:
                scope = json.load(fh).get("baf_scope") or {}
        except (OSError, ValueError):
            scope = {}

    # median on-target log2 per window from the 1 kb bins
    log2_by_gene = defaultdict(list)
    bins_by_gene = defaultdict(list)
    for b in bins:
        v = _float(b.get("log2_ratio"))
        if v is not None:
            log2_by_gene[b.get("gene")].append(v)
            bins_by_gene[b.get("gene")].append({"start": _int(b.get("start")) or 0,
                                                "end": _int(b.get("end")) or 0, "log2": v})
    sites_by_gene = defaultdict(list)
    for s in _read_tsv(os.path.join(d, "%s.genebaf.sites.tsv" % sample)):
        af = _float(s.get("alt_fraction"))
        pos = _int(s.get("pos"))
        if af is not None and pos is not None:
            sites_by_gene[s.get("gene")].append({"pos": pos, "af": af})

    cn_by_gene = {}
    for g in (ichorkaryo_genes or []):
        cn_by_gene[g.get("gene")] = g

    by_chrom = defaultdict(lambda: {"windows": [], "arms": []})
    for w in win:
        gene = w.get("gene")
        cnrec = cn_by_gene.get(gene, {})
        state = (w.get("allelic_state") or "").strip()
        obs = _float(w.get("observed_het_fraction"))
        exp = _float(w.get("expected_het_fraction"))
        lg = log2_by_gene.get(gene)
        row = {
            "gene": gene, "chrom": w.get("chrom"),
            "start": _int(w.get("start")), "end": _int(w.get("end")),
            "size_kb": ((_int(w.get("end")) or 0) - (_int(w.get("start")) or 0)) / 1000.0,
            "band": cnrec.get("band", ""),
            "mean_depth": _float(w.get("mean_depth")),
            "log2_ontarget": round(statistics.median(lg), 3) if lg else None,
            "cn": cnrec.get("cn"), "cn_event": cnrec.get("event", ""),
            "cell_fraction": cnrec.get("cell_fraction"),
            "baf_scope": w.get("baf_scope", ""),
            "allelic_state": state,
            "n_het": _int(w.get("n_het")), "n_usable": _int(w.get("n_usable")),
            "median_dev": _float(w.get("median_dev_from_half")),
            "observed_het": obs, "expected_het": exp,
        }
        row["reading"] = _reading(row["cn_event"], state, obs, exp)
        row["genes_in_window"] = _genes_in_window(
            gene_model, row["chrom"], row["start"] or 0, row["end"] or 0)
        row["svg"] = (_window_svg(gene, row["start"] or 0, row["end"] or 0,
                                  bins_by_gene.get(gene, []), sites_by_gene.get(gene, []),
                                  row["genes_in_window"])
                      if row["baf_scope"] == "assessed" else "")
        by_chrom[row["chrom"]]["windows"].append(row)

    for a in arms:
        row = {
            "arm": a.get("arm"), "verdict": a.get("verdict", ""),
            "n_windows": _int(a.get("n_windows")), "genes": a.get("genes", ""),
            "n_het": _int(a.get("n_het")), "n_usable": _int(a.get("n_usable")),
            "ratio": _float(a.get("ratio")), "p_value": _float(a.get("p_value")),
            "span_start": _int(a.get("span_start")), "span_end": _int(a.get("span_end")),
            "covered_fraction": _float(a.get("covered_fraction")),
            "span_fraction": _float(a.get("span_fraction")),
            "arm_bp": _int(a.get("arm_bp")),
        }
        by_chrom[a.get("chrom")]["arms"].append(row)

    for chrom in by_chrom:
        by_chrom[chrom]["windows"].sort(key=lambda r: r["start"] or 0)
    result.update({
        "found": True,
        "chromosomes": {c: by_chrom[c] for c in CHROM_ORDER if c in by_chrom},
        "scope": scope,
        "n_windows": len(win),
        "n_assessed": sum(1 for w in win if (w.get("baf_scope") or "") == "assessed"),
        "n_imbalanced": sum(1 for w in win if (w.get("allelic_state") or "") == "imbalanced"),
    })
    return result
