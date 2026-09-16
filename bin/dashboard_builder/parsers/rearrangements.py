"""
rearrangements.py - dashboard parser for the merged, graded rearrangement table.

Reads <effective_dir>/translocations/<sample>.translocations.tsv, the output
of merge_translocations.py (one row per junction after positional collapse),
and organises it for the Rearrangements tab in three levels:

    event group   - every junction that touches the same anchor locus
                    (MYC, IGH, IGK, IGL/IGLL5). MYC in particular is often
                    rearranged with several partners at once.
    partner locus - within a group, junctions whose two breakends both fall
                    within `locus_tol` bp of each other are one partner
                    locus. This is the row the reader sees.
    breakend      - the individual junctions, listed under their locus.

A junction with no anchor is its own group and its own locus.

Two tables come out of this: `reportable` (any member has reportable=yes)
and `other` (everything else, including flagged junctions). Each keeps the
same structure so the reader learns one layout.

The tolerance is a display parameter, not a biological claim. Partner-side
breakpoints of one MYC rearrangement cluster within the TAD that holds the
partner's super-enhancer, and the MYC-side hotspots lie within a 2.5 Mb
span, so 1 Mb on each side collects one event without reaching the next
locus on the same chromosome.

Also reads <sample>.mm_annotated.tsv only to count and list the non-
rearrangement structural variants (DEL/INS/INV/DUP) that carry a tier, so a
defining intrachromosomal event is never hidden.

Standard library only.
"""

import csv
import json
import os
import re
from collections import defaultdict

ANCHOR_PRIORITY = ["MYC", "IGH", "IGK", "IGL"]
IG_TOKENS = {"IGH", "IGK", "IGL", "IGLL5"}
TRANSLOCATION_TYPES = {"TRA", "BND"}
FILTER_RANK = {"PASS": 0, "": 0}

# Plain-language reading key. Rendered once above the tables.
READING_KEY = [
    ("Reportable",
     "A junction is reportable if it touches IGH, IGK, IGL or MYC, or if the "
     "pair is a named myeloma translocation. A flagged junction is not "
     "reportable. Everything else is of unknown significance and is listed "
     "separately."),
    ("Tier",
     "defining: a primary Ig translocation that defines the disease entity "
     "(WHO 5th). secondary: MYC and other progression events. rare: published "
     "partners seen in under 1% of cases. A trailing ? means only one side "
     "was resolved to a gene."),
    ("Flags",
     "A reason to doubt a junction; the row is kept but not reported. "
     "MYC with a non-Ig partner: in myeloma, MYC is activated by being placed "
     "next to an immunoglobulin enhancer, so a MYC partner that is not IGH, "
     "IGK or IGL is unusual. IGH V region with an off-panel partner: the V "
     "region is a stretch of about 120 near-identical gene segments; reads "
     "from one segment often align to another, which produces false "
     "breakends, so a breakend there whose partner is not a known gene is "
     "most likely mapping error. Multi-partner breakend: the same position "
     "joins three or more different places in the genome; a real breakend "
     "joins one."),
    ("Ig region",
     "Where in the locus a breakend sits: C_switch (constant/switch region, "
     "where the primary translocations arise), J, D, V (see above). "
     "Informational only."),
    ("Reads and callers",
     "Read counts come from each caller's own output. Worst filter is the "
     "least favourable verdict among the callers that saw the junction: PASS; "
     "q5 = CuteSV's quality score below 5, which every two- or three-read "
     "call receives, so read it as low support rather than error; GT = "
     "Sniffles could not assign a genotype, likewise a low-support signal; "
     "PREDICTED_NOISE = SAVANA's classifier scored the call as noise. Use "
     "these to rank, not to exclude: true translocations at low tumour "
     "fraction carry these labels."),
    ("nanomonsv",
     "An independent caller run with a panel of normal genomes, used to "
     "confirm rather than to vote. confirmed: it found the same junction. "
     "not reported: it did not; at low support this means nothing either "
     "way. A filter name: it saw the junction and rejected it for that "
     "reason."),
    ("Grouping",
     "Junctions sharing an anchor (typically MYC, which is often rearranged "
     "with several partners at once) form one event group. Within a group, "
     "breakends whose two sides fall within 1 Mb of each other on both "
     "chromosomes are one partner locus. Open a partner locus to see its "
     "breakends."),
]

COLUMN_KEY = [
    ("Anchor / partner", "the anchor locus and the partner locus of the junction; coordinates on T2T-CHM13v2.0 with the nearest gene or panel window, cytoband and Ig sub-region"),
    ("Named pair / tier", "dictionary name and tier; entity per WHO 5th where the pair defines one"),
    ("Reads", "maximum supporting reads across callers and across the breakends of the locus"),
    ("Callers", "ensemble callers (Sniffles, CuteSV, Severus, SAVANA) that reported any breakend of the locus"),
    ("Worst filter", "least favourable per-caller FILTER on the best-supported breakend; see key"),
    ("nanomonsv", "confirmation status and nanomonsv's own read count"),
    ("Flag", "hover for the reason; a flagged row is never reportable"),
    ("Breakends", "number of distinct junctions collapsed into this locus; open the row to see them"),
]


def _find(effective_dir, sample, suffix):
    direct = os.path.join(effective_dir, "translocations", "%s.%s" % (sample, suffix))
    if os.path.isfile(direct):
        return direct
    for root, _d, files in os.walk(effective_dir):
        for name in files:
            if name == "%s.%s" % (sample, suffix):
                return os.path.join(root, name)
    return None


def _load_igv_manifest(effective_dir, sample):
    """event_id -> {"a": relpath, "b": relpath} from igv_snapshots.py."""
    manifest_path = None
    for root, _d, files in os.walk(effective_dir):
        for name in files:
            if name == "%s.translocations.manifest.json" % sample:
                manifest_path = os.path.join(root, name)
                break
        if manifest_path:
            break
    if not manifest_path:
        return {}
    try:
        with open(manifest_path) as fh:
            manifest = json.load(fh)
    except (OSError, ValueError):
        return {}
    base = os.path.dirname(manifest_path)
    lookup = {}
    for event in manifest.get("events", []):
        eid = event.get("event_id")
        if not eid:
            continue
        entry = {}
        for side in ("a", "b"):
            html = (event.get(side) or {}).get("html")
            if not html:
                continue
            absolute = html if os.path.isabs(html) else os.path.join(base, html)
            if os.path.isfile(absolute):
                entry[side] = os.path.relpath(absolute, effective_dir)
        if entry:
            lookup[eid] = entry
    return lookup


def _int(value):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _tokens(label):
    return {t.strip().upper() for t in re.split(r"[/,+ ]", label or "") if t.strip()}


def _anchor_of(row):
    """Which anchor this junction belongs to and which side carries it.

    Returns (anchor_label, side) with side in {"a", "b"} or (None, None).
    MYC first: a MYC::IGH junction is a MYC event (secondary), not an IGH one.
    """
    ga, gb = _tokens(row.get("gene_a")), _tokens(row.get("gene_b"))
    for anchor in ANCHOR_PRIORITY:
        if anchor in ga or (anchor == "IGL" and "IGLL5" in ga):
            return anchor, "a"
        if anchor in gb or (anchor == "IGL" and "IGLL5" in gb):
            return anchor, "b"
    return None, None


def _filter_rank(value):
    v = (value or "").strip()
    if v in FILTER_RANK:
        return 0
    return 1 + len(v.split(";"))


def _mb(pos):
    return "%.2f" % (pos / 1e6)


def _span_label(chrom, lo, hi):
    if lo == hi:
        return "%s:%s" % (chrom, "{:,}".format(lo))
    return "%s:%s–%s Mb" % (chrom, _mb(lo), _mb(hi))


def _cluster(members, tol):
    """Single-linkage on both sides within tol bp. members carry
    'apos'/'ppos' (anchor side / partner side, or a/b when unanchored)."""
    members = sorted(members, key=lambda m: (m["pchrom"], m["ppos"], m["apos"]))
    clusters = []
    for m in members:
        placed = None
        for c in clusters:
            if c["pchrom"] != m["pchrom"]:
                continue
            if any(abs(m["ppos"] - x["ppos"]) <= tol and abs(m["apos"] - x["apos"]) <= tol
                   for x in c["members"]):
                placed = c
                break
        if placed is None:
            clusters.append({"pchrom": m["pchrom"], "members": [m]})
        else:
            placed["members"].append(m)
    return clusters


def _best(members, key, default=0):
    vals = [m[key] for m in members if m[key] is not None]
    return max(vals) if vals else default


def _mode(values):
    counts = defaultdict(int)
    for v in values:
        if v:
            counts[v] += 1
    if not counts:
        return ""
    return max(counts.items(), key=lambda kv: (kv[1], -len(kv[0])))[0]


def parse(effective_dir, sample, locus_tol=1_000_000, igv_flank=5000):
    def not_found(reason):
        return {"found": False, "reason": reason, "searched": str(effective_dir),
                "reportable": [], "other": [], "n_junctions": 0, "summary": {},
                "reading_key": READING_KEY, "column_key": COLUMN_KEY,
                "other_svs": [], "locus_tol": locus_tol}

    path = _find(effective_dir, sample, "translocations.tsv")
    if not path:
        return not_found("no %s.translocations.tsv under this directory" % sample)

    with open(path, newline="") as fh:
        rows = [{k: (v or "").strip() for k, v in r.items()}
                for r in csv.DictReader(fh, delimiter="\t")]

    igv = _load_igv_manifest(effective_dir, sample)

    # ---- per-junction records ------------------------------------------
    junctions = []
    for r in rows:
        anchor, side = _anchor_of(r)
        a, b = ("a", "b") if side in (None, "a") else ("b", "a")
        pos_a, pos_b = _int(r.get("pos_" + a)), _int(r.get("pos_" + b))
        if pos_a is None or pos_b is None:
            continue
        ids = [r.get("sv_id")] + [x for x in (r.get("merged_sv_ids") or "").split(",") if x]
        pages = {}
        for eid in ids:
            if eid in igv:
                pages = igv[eid]
                break
        detail_keys = ("sv_id", "merged_sv_ids", "n_merged", "support_representative",
                       "support_sniffles", "support_cutesv", "support_severus", "support_savana",
                       "filter_sniffles", "filter_cutesv", "filter_severus", "filter_savana",
                       "support_nanomonsv", "filter_nanomonsv", "supp_vec", "pon_freq",
                       "known_mm_pair", "entity", "tier", "known_freq", "match_quality",
                       "anchor", "anchor_class", "flag", "dict_notes",
                       "ig_region_a", "ig_region_b", "gene_a_dist", "gene_b_dist", "filter")
        junctions.append({
            "detail": [(k, r.get(k)) for k in detail_keys if (r.get(k) or "").strip()],
            "sv_id": r.get("sv_id"),
            "anchor": anchor,
            "achrom": r.get("chrom_" + a), "apos": pos_a,
            "agene": r.get("gene_" + a), "aband": r.get("band_" + a),
            "aig": r.get("ig_region_" + a),
            "pchrom": r.get("chrom_" + b), "ppos": pos_b,
            "pgene": r.get("gene_" + b), "pband": r.get("band_" + b),
            "pig": r.get("ig_region_" + b),
            "sv_type": (r.get("sv_type") or "").upper(),
            "reads": _int(r.get("support_reads")),
            "callers": [c.strip() for c in (r.get("callers") or "").split(",") if c.strip()],
            "worst_filter": r.get("filter_worst") or "",
            "nano_conf": r.get("nanomonsv_confirmed") or "",
            "nano_filter": r.get("filter_nanomonsv") or "",
            "nano_reads": _int(r.get("support_nanomonsv")),
            "known_pair": r.get("known_mm_pair") or "",
            "tier": r.get("tier") or "",
            "entity": r.get("entity") or "",
            "known_freq": r.get("known_freq") or "",
            "flag": r.get("flag") or "",
            "reportable": (r.get("reportable") or "").lower() == "yes",
            "n_merged": _int(r.get("n_merged")) or 1,
            "igv": {
                "html_a": pages.get("a"), "html_b": pages.get("b"),
                "point_a": "%s:%s" % (r.get("chrom_a"), r.get("pos_a")),
                "point_b": "%s:%s" % (r.get("chrom_b"), r.get("pos_b")),
                "label_a": r.get("gene_a") or r.get("chrom_a"),
                "label_b": r.get("gene_b") or r.get("chrom_b"),
                "event_id": r.get("sv_id"),
            },
        })

    # ---- group -> locus ---------------------------------------------------
    by_anchor = defaultdict(list)
    for j in junctions:
        key = j["anchor"] or ("pair:%s|%s" % (j["agene"], j["pgene"]))
        by_anchor[key].append(j)

    def build_locus(cluster, anchor):
        ms = sorted(cluster["members"], key=lambda m: (-(m["reads"] or 0), m["ppos"]))
        best = ms[0]
        apos = [m["apos"] for m in ms]
        ppos = [m["ppos"] for m in ms]
        callers = sorted({c for m in ms for c in m["callers"]},
                         key=lambda c: ["Sniffles", "CuteSV", "Severus", "SAVANA"].index(c)
                         if c in ("Sniffles", "CuteSV", "Severus", "SAVANA") else 9)
        nano_vals = [m["nano_conf"] for m in ms]
        nano = "yes" if "yes" in nano_vals else ("no" if "no" in nano_vals else "")
        nano_filter = _mode([m["nano_filter"] for m in ms if m["nano_conf"] == "no"])
        tiered = [m for m in ms if m["tier"]]
        graded = tiered[0] if tiered else None
        flagged = [m for m in ms if m["flag"]]
        return {
            "anchor": anchor,
            "anchor_label": _span_label(best["achrom"], min(apos), max(apos)),
            "anchor_gene": _mode([m["agene"] for m in ms]),
            "anchor_band": _mode([m["aband"] for m in ms]),
            "anchor_ig": _mode([m["aig"] for m in ms]),
            "partner_label": _span_label(best["pchrom"], min(ppos), max(ppos)),
            "partner_gene": _mode([m["pgene"] for m in ms]),
            "partner_band": _mode([m["pband"] for m in ms]),
            "partner_ig": _mode([m["pig"] for m in ms]),
            "partner_chrom": best["pchrom"],
            "n_breakends": len(ms),
            "n_records": sum(m["n_merged"] for m in ms),
            "reads": _best(ms, "reads"),
            "callers": callers,
            "n_callers": len(callers),
            "worst_filter": best["worst_filter"],
            "nano": nano,
            "nano_filter": nano_filter,
            "nano_reads": _best(ms, "nano_reads", None),
            "known_pair": graded["known_pair"] if graded else "",
            "tier": graded["tier"] if graded else "",
            "entity": graded["entity"] if graded else "",
            "known_freq": graded["known_freq"] if graded else "",
            "flag": flagged[0]["flag"] if (flagged and len(flagged) == len(ms)) else "",
            "reportable": any(m["reportable"] for m in ms),
            "igv": best["igv"],
            "detail": best["detail"],
            "members": ms,
            "sort_key": (0 if any(m["reportable"] for m in ms) else 1,
                         0 if graded else 1, -(_best(ms, "reads") or 0)),
        }

    groups = []
    for key, js in by_anchor.items():
        anchor = js[0]["anchor"]
        loci = [build_locus(c, anchor) for c in _cluster(js, locus_tol)]
        loci.sort(key=lambda l: l["sort_key"])
        groups.append({
            "anchor": anchor,
            "anchor_display": {"MYC": "MYC", "IGH": "IGH", "IGK": "IGK", "IGL": "IGL/IGLL5"}.get(anchor, ""),
            "loci": loci,
            "n_loci": len(loci),
            "reportable": any(l["reportable"] for l in loci),
            "sort_key": min(l["sort_key"] for l in loci),
        })

    def split(pred):
        out = []
        for g in groups:
            loci = [l for l in g["loci"] if pred(l)]
            if loci:
                gg = dict(g)
                gg["loci"] = loci
                gg["n_loci"] = len(loci)
                gg["sort_key"] = min(l["sort_key"] for l in loci)
                out.append(gg)
        out.sort(key=lambda g: g["sort_key"])
        return out

    reportable = split(lambda l: l["reportable"])
    other = split(lambda l: not l["reportable"])

    # ---- non-rearrangement SVs carrying a tier, from mm_annotated -------
    other_svs = []
    n_all = 0
    ann = _find(effective_dir, sample, "mm_annotated.tsv")
    if ann:
        with open(ann, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                n_all += 1
                svt = (r.get("sv_type") or "").upper()
                if svt in TRANSLOCATION_TYPES:
                    continue
                if r.get("chrom_a") != r.get("chrom_b"):
                    continue
                if not (r.get("tier") or "").strip():
                    continue
                other_svs.append({
                    "sv_id": r.get("sv_id"), "sv_type": svt,
                    "locus": "%s:%s-%s" % (r.get("chrom_a"), r.get("pos_a"), r.get("pos_b")),
                    "genes": "%s / %s" % (r.get("gene_a"), r.get("gene_b")),
                    "reads": r.get("support_reads"), "callers": r.get("callers"),
                    "tier": r.get("tier"), "known_pair": r.get("known_mm_pair"),
                    "filter": r.get("filter_worst"),
                })

    all_loci = [l for g in groups for l in g["loci"]]
    defining = [l for l in all_loci if l["tier"].startswith("defining")]
    summary = {
        "n_junctions": len(junctions),
        "n_all_svs": n_all,
        "n_loci": len(all_loci),
        "n_reportable": sum(1 for l in all_loci if l["reportable"]),
        "n_defining": len(defining),
        "defining_names": sorted({l["known_pair"] for l in defining if l["known_pair"]}),
        "n_confirmed": sum(1 for l in all_loci if l["nano"] == "yes"),
        "n_flagged": sum(1 for l in all_loci if l["flag"]),
        "n_groups_multi": sum(1 for g in groups if g["n_loci"] > 1 and g["anchor"]),
        "has_nanomonsv": any(m["nano_conf"] for m in junctions),
    }

    return {
        "found": True,
        "path": path, "filename": os.path.basename(path),
        "reportable": reportable, "other": other,
        "other_svs": other_svs,
        "summary": summary,
        "reading_key": READING_KEY, "column_key": COLUMN_KEY,
        "locus_tol": locus_tol,
        "has_igv": bool(igv),
    }
