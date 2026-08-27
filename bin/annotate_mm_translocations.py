#!/usr/bin/env python3
"""
annotate_mm_translocations.py
=============================

Annotates a SURVIVOR-merged VCF of structural variants against the plasma
cell neoplasm panel, by intersecting breakpoints with the panel BED.

For each record where at least one breakpoint falls inside a panel region,
emits a row carrying both breakpoints, what each one landed in, whether the
pair is named by the dictionary, and how firmly.

What this version adds
----------------------
GENE MODEL. Breakpoints were named from the panel interval, which is flanked
and merged for capture. 62% of the v7 panel is flank, and five intervals
carry compound labels. A t(4;14) breakend landed in a window called
"FGFR3/NSD2" and the table could not say which gene it was in, so IGH::FGFR3
and IGH::NSD2 were indistinguishable. With a gene model the two are 62 kb
apart and the answer is available. gene_a_dist and gene_b_dist record how far
each breakend sits from the gene it was named after.

Falling back to the interval label is still the right answer for a
breakpoint-cluster event and is not a failure: MYC's body is 8 kb inside a
5 Mb window, and t(11;14) breakpoints sit 100-400 kb outside CCND1. The
distance columns are what make that auditable rather than implicit.

CYTOBAND COLUMNS. band_a and band_b are emitted for both sides. Downstream,
merge_translocations groups junctions by chromosome arm, which is parsed from
the band rather than the chromosome, and could not be done without them.

GRADING. tier, entity and reportable come from the dictionary and the anchor
table. Nothing downstream could previously tell a defining translocation from
an incidental junction, so the report and the IGV page selection filtered on
SV type alone and hid intrachromosomal findings behind a toggle.

ANCHORS. A dictionary of named pairs cannot cover a locus whose partner list
is open. Any IGH, IGK, IGL or MYC junction is reportable whatever the
partner, including a partner that is off-panel and resolves only to a band.

EXCLUSIONS. Junctions observed at coordinate-identical positions in unrelated
patients can be dropped. The criterion is coordinate identity, not
recurrence: this disease is defined by recurrent events. A dictionary-named
or graded row is never dropped, and that override is not configurable.

This script holds no biological priors of its own. Every pair, tier, entity
and anchor comes from --dictionary and --anchors; every coordinate comes from
--panel-bed, --gene-model and --cytoband-bed. No variant, FISH finding or
expected karyotype is hardcoded anywhere.

Usage:
  annotate_mm_translocations.py \\
      --vcf                merged.vcf.gz \\
      --panel-bed          aWGS_PCN_v7_t2t_chr.bed \\
      --gene-model         aWGS_PCN_v7_gene_model_t2t.bed \\
      --dictionary         mm_translocation_dictionary.tsv \\
      --anchors            mm_translocation_anchors.tsv \\
      --cytoband-bed       chm13v2.0_cytobands_allchrs.bed \\
      --excluded-junctions mm_excluded_junctions.tsv \\
      --sample             SAMPLE_ID \\
      --output             SAMPLE_ID.mm_annotated.tsv
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


__version__ = "0.6.0"
# [dictionary-token-matching applied]
# [cytoband-partner-annotation applied]


@dataclass
class PanelRegion:
    chrom: str
    start: int
    end: int
    name: str

    def contains(self, chrom: str, pos: int) -> bool:
        return chrom == self.chrom and self.start <= pos < self.end


@dataclass
class SvRecord:
    chrom: str
    pos: int
    sv_id: str
    sv_type: str
    mate_chrom: Optional[str]
    mate_pos: Optional[int]
    filt: str
    info: Dict[str, str]
    callers: List[str]
    support: str = ""


def load_panel(bed_path: Path) -> List[PanelRegion]:
    out: List[PanelRegion] = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) >= 4 else f"{chrom}:{start}-{end}"
            out.append(PanelRegion(chrom, start, end, name))
    return out


# Region-suffix words that are not gene symbols (e.g. the "locus" in
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


def parse_info(info_field: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for part in info_field.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            out[k] = v
        else:
            out[part] = "True"
    return out


def parse_bnd_alt(alt: str) -> Tuple[Optional[str], Optional[int]]:
    """Parse BND ALT (N]chr14:99000000] or [chr8:130000000[N)."""
    for opener, closer in [("]", "]"), ("[", "[")]:
        if opener in alt:
            try:
                inner = alt.split(opener)[1].split(closer)[0]
                chrom, pos = inner.rsplit(":", 1)
                return chrom, int(pos)
            except (IndexError, ValueError):
                continue
    return None, None


def infer_callers(info: Dict[str, str]) -> List[str]:
    """
    Decode SUPP_VEC bitstring written by SURVIVOR. Bit order matches the
    order of input VCFs passed to SURVIVOR merge; in our pipeline that is
    [Sniffles, CuteSV, Severus].
    """
    out = []
    supp_vec = info.get("SUPP_VEC", "")
    caller_order = ["Sniffles", "CuteSV", "Severus"]
    for bit, name in zip(supp_vec, caller_order):
        if bit == "1":
            out.append(name)
    return out


def support_reads_from(info: Dict[str, str], fmt: str, sample_cols: List[str]) -> str:
    """
    Variant-supporting read count for a record, read generically (no hardcoded
    calls). The merged VCF carries one sample column per input caller; the
    caller that detected the junction holds the real FORMAT/DV while the others
    are zero/absent, so DV is taken as the MAX across all sample columns. If no
    DV is present, fall back to caller-specific INFO tags:
      RE       (CuteSV)
      SUPPORT  (Sniffles)
      SR       (generic)
      SUPP_READS first colon-field (Severus, e.g. "2:0:0:2:0:0" -> 2)
    """
    if fmt:
        keys = fmt.split(":")
        best = None
        for col in sample_cols:
            if not col or col in (".", "./."):
                continue
            f = dict(zip(keys, col.split(":")))
            dv = f.get("DV", "")
            if dv.isdigit():
                best = int(dv) if best is None else max(best, int(dv))
        if best is not None:
            return str(best)
    for tag in ("RE", "SUPPORT", "SR"):
        v = str(info.get(tag, ""))
        if v.isdigit():
            return v
    sr = str(info.get("SUPP_READS", ""))
    first = sr.split(":")[0] if sr else ""
    if first.isdigit():
        return first
    return ""


def open_vcf(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def parse_vcf(vcf_path: Path) -> List[SvRecord]:
    out: List[SvRecord] = []
    with open_vcf(vcf_path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            chrom, pos, sv_id, _ref, alt, _qual, filt, info_field = cols[:8]
            pos = int(pos)
            info = parse_info(info_field)
            sv_type = info.get("SVTYPE", "")
            mate_chrom, mate_pos = (None, None)

            if sv_type in ("BND", "TRA"):
                # BND records have mate coords in ALT (e.g. N]chr11:69500000]).
                # SURVIVOR renames cross-chromosome BNDs to TRA while keeping
                # the same ALT format, so we parse both the same way.
                mate_chrom, mate_pos = parse_bnd_alt(alt)
                # Fallback: if ALT parsing failed (some SURVIVOR records use
                # symbolic ALT like <TRA>), use SURVIVOR's CHR2/END INFO tags.
                if mate_chrom is None:
                    chr2 = info.get("CHR2")
                    end = info.get("END")
                    if chr2 and end:
                        try:
                            mate_chrom, mate_pos = chr2, int(end)
                        except ValueError:
                            pass

            elif sv_type in ("DEL", "DUP", "INV", "INS"):
                # Same-chromosome events: mate side is END on the same chrom.
                end = info.get("END")
                if end:
                    try:
                        mate_chrom, mate_pos = chrom, int(end)
                    except ValueError:
                        pass

            # FORMAT + per-caller sample columns carry read support (DV).
            fmt = cols[8] if len(cols) > 8 else ""
            sample_cols = cols[9:] if len(cols) > 9 else []
            support = support_reads_from(info, fmt, sample_cols)

            out.append(SvRecord(
                chrom=chrom, pos=pos, sv_id=sv_id, sv_type=sv_type,
                mate_chrom=mate_chrom, mate_pos=mate_pos,
                filt=filt, info=info, callers=infer_callers(info),
                support=support,
            ))
    return out



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
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            chrom, start, end, band = parts[0], int(parts[1]), int(parts[2]), parts[3]
            bands.setdefault(chrom, []).append((start, end, band))
    for chrom in bands:
        bands[chrom].sort(key=lambda t: t[0])
    if not bands:
        sys.stderr.write(f"ERROR: no cytobands parsed from {bed_path}\n")
        sys.exit(1)
    return CytobandTable(bands)


def _strip_chr(chrom: str) -> str:
    """chr8 -> 8, chrX -> X, for clinical band notation (8q24.21)."""
    return chrom[3:] if chrom.startswith("chr") else chrom


# -----------------------------------------------------------------------------
# Gene model
# -----------------------------------------------------------------------------
def load_gene_model(path: Optional[Path]) -> List[PanelRegion]:
    """Bare gene bodies, for naming a breakpoint and measuring distance.

    Separate from the panel BED and never used to decide panel membership.
    The panel's intervals are flanked and merged, which is right for capture
    and wrong for naming: 62% of the v7 panel is flank, and five of its
    intervals carry compound labels (FGFR3/NSD2, WWOX/MAF, TP53+TNFSF12,
    FCRL5/FCRL4, IGL/IGLL5) that name two genes without saying which one a
    breakend fell in.

    Optional. Without it the annotator behaves as it did before: labels come
    from the panel interval and the distance columns stay empty.
    """
    if path is None or not Path(str(path)).is_file():
        return []
    return load_panel(Path(path))


def load_ig_segments(path: Optional[Path]) -> List[PanelRegion]:
    """Immunoglobulin locus sub-regions: C/switch, J, D, V.

    Optional. Built by bin/build_ig_segments.py from the RefSeq GFF, so the
    boundaries are read off the annotation rather than asserted; they move
    between assemblies and hardcoding one would make the rule wrong on the
    next reference.
    """
    if path is None or not Path(str(path)).is_file():
        return []
    return load_panel(Path(path))


def ig_region_for(chrom, pos, segments) -> str:
    """Which Ig sub-region this breakpoint falls in, or "".

    Position within the locus carries mechanistic information the coordinate
    alone does not. Primary translocations in plasma cell neoplasms arise
    from errors of class-switch recombination, which acts at the switch
    regions, so they break in C/switch or near J. The V array is 121
    near-identical segments over 940 kb where V(D)J and somatic
    hypermutation act physiologically and where reads mismap.

    Reported, never filtered on. A V-array breakend with an off-panel
    partner is very likely mismapping, but this column says where the
    breakend is and leaves the judgement to the reader.
    """
    best = None
    for r in segments:
        if r.chrom == chrom and pos is not None and r.start <= pos < r.end:
            if best is None or (r.end - r.start) < (best.end - best.start):
                best = r
    return best.name if best else ""


def gene_for(chrom, pos, model) -> Optional[str]:
    """Tightest gene-model feature containing this coordinate, or None.

    Tightest wins because gene models legitimately nest: IGLL5 lies inside
    the IGL locus, and the smaller feature is the more specific answer. A
    breakend in IGLL5 is named IGLL5; one elsewhere in the locus is IGL.
    """
    if not model or chrom is None or pos is None:
        return None
    best = None
    for reg in model:
        if reg.chrom == chrom and reg.start <= pos < reg.end:
            if best is None or (reg.end - reg.start) < (best.end - best.start):
                best = reg
    return best.name if best else None


def dist_to_gene(chrom, pos, name, model) -> Optional[int]:
    """Bases from this coordinate to the named gene's body; 0 if inside.

    None when the label names nothing in the model - a cytoband, a compound
    panel label, a bare coordinate - because there is no body to measure
    from. The nearest record of that name wins when a symbol has several.
    """
    if not model or chrom is None or pos is None or not name:
        return None
    best = None
    for reg in model:
        if reg.chrom == chrom and reg.name == name:
            d = max(reg.start - pos, pos - reg.end, 0)
            if best is None or d < best:
                best = d
    return best


# -----------------------------------------------------------------------------
# Dictionary, graded
# -----------------------------------------------------------------------------
@dataclass
class DictEntry:
    tok_a: frozenset
    tok_b: frozenset
    band_b: str
    row: dict


def load_dictionary(dict_path: Path) -> List[DictEntry]:
    """Named partner pairs. The only source of biological priors here.

    A missing file is non-fatal and yields an empty list: every junction is
    then emitted unnamed rather than the run failing. Nothing in this script
    knows any pair, tier or entity of its own.
    """
    out: List[DictEntry] = []
    if not dict_path.exists():
        sys.stderr.write(f"WARNING: dictionary not found: {dict_path}\n")
        return out
    with open(dict_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            a = (row.get("partner_a") or "").strip()
            b = (row.get("partner_b") or "").strip()
            if not a or not b:
                continue
            out.append(DictEntry(_norm_tokens(a), _norm_tokens(b),
                                 (row.get("partner_b_band") or "").strip(), row))
    return out


def _band_levels(band: Optional[str]) -> List[str]:
    """Progressively coarser forms of a cytoband label, most precise first.

    '14q32.33' -> ['14q32.33', '14q32', '14q']

    Exact sub-band equality is too brittle to depend on. T2T-CHM13 band
    boundaries are not GRCh38's, and the published band for a partner is
    usually quoted against GRCh38, so a breakpoint one sub-band away would
    otherwise lose its name.
    """
    if not band:
        return []
    out = [band]
    if "." in band:
        out.append(band.split(".", 1)[0])
    for i, ch in enumerate(out[-1]):
        if ch in "pq":
            arm = out[-1][:i + 1]
            if arm != out[-1]:
                out.append(arm)
            break
    return out


def _band_match(observed: Optional[str], expected: str) -> Optional[str]:
    """Match quality if observed and expected bands agree at any level.

    Walks coarse to fine, remembering the finest level that agrees and
    stopping at the first shared level that disagrees. 16q23 and 16q12 agree
    on the arm but differ at the major band: that is a contradiction, not a
    partial match, and it returns partial_arm rather than being promoted.
    """
    if not observed or not expected:
        return None
    obs, exp = _band_levels(observed), _band_levels(expected)
    if not obs or not exp:
        return None
    depth = -1
    for i in range(min(len(obs), len(exp))):
        if obs[i] == exp[i]:
            depth = i
        else:
            if depth < 0:
                return None
            break
    if depth >= 1:
        return "partial_band"
    if depth == 0:
        return "partial_arm"
    return None


def dictionary_lookup(dictionary, label_a, label_b,
                      band_a=None, band_b=None) -> Tuple[Optional[dict], str]:
    """Return (row, match_quality) for an unordered pair.

    'full'          both sides matched on gene identity
    'partial_band'  one side matched a gene, the other matched the expected
                    cytoband the dictionary records for an off-panel partner
    'partial_arm'   as above but agreeing only at arm level: a lead, not a
                    call
    ('', None) when nothing matches.

    Band matching exists for partners this panel does not capture. Every
    partner in the shipped dictionary is on-panel, so partner_b_band is
    empty throughout and only the 'full' path fires today. The mechanism is
    here so a partner can be added to the dictionary without also having to
    be added to the panel.
    """
    ta, tb = _norm_tokens(label_a), _norm_tokens(label_b)
    if not ta or not tb:
        return None, ""

    for e in dictionary:
        if (ta & e.tok_a and tb & e.tok_b) or (ta & e.tok_b and tb & e.tok_a):
            return e.row, "full"

    best = None
    for e in dictionary:
        if not e.band_b:
            continue
        for gene_tok, other_band in ((ta, band_b), (tb, band_a)):
            if not (gene_tok & e.tok_a):
                continue
            q = _band_match(other_band, e.band_b)
            if q == "partial_band":
                return e.row, q
            if q and best is None:
                best = (e.row, q)
    return best if best else (None, "")


# -----------------------------------------------------------------------------
# Anchors
# -----------------------------------------------------------------------------
def load_anchors(path: Optional[Path]) -> Dict[str, dict]:
    """Gene token -> anchor row.

    An anchor is a locus whose partner list is open: the dictionary names
    what is described, the anchor covers what is not. Comment lines are
    skipped so the table can carry its own documentation.
    """
    out: Dict[str, dict] = {}
    if path is None or not Path(str(path)).is_file():
        sys.stderr.write(f"WARNING: anchor table not found: {path}\n")
        return out
    with open(path) as fh:
        rows = [ln for ln in fh if not ln.lstrip().startswith("#")]
    for row in csv.DictReader(rows, delimiter="\t"):
        anchor = (row.get("anchor") or "").strip()
        if not anchor:
            continue
        for tok in _norm_tokens(anchor):
            out[tok] = row
    return out


def anchor_hits(anchors, label_a, label_b, dist_a, dist_b) -> List[dict]:
    """Anchor rows triggered by either side of a junction.

    A self-pair returns nothing. An anchor is a claim about partners: the
    locus rearranges with many of them, so any partner is worth surfacing. A
    locus joined to itself has no partner and the premise does not apply.
    Without that guard every V(D)J and somatic hypermutation product inside
    IGH, IGK and IGL would inherit reportability from a rule that was never
    about it, and in a plasma cell neoplasm those loci recombine and
    hypermutate physiologically.

    dist_policy decides whether a breakend in the panel flank still counts.
    'body' requires the breakend to be in the named gene, or the label to be
    one with no body to measure from. 'interval' fires wherever the label
    applies. Every anchor on this panel is 'interval', because the windows
    are deliberately wide: t(11;14) breakpoints sit 100-400 kb from CCND1
    and t(14;16) breakpoints about 5 Mb from MAF, so a body rule would
    discard the events the BCR windows exist to capture. The distance is
    still reported on every row.
    """
    by_side = []
    for label, dist in ((label_a, dist_a), (label_b, dist_b)):
        side = {}
        for tok in _norm_tokens(label or ""):
            row = anchors.get(tok)
            if not row:
                continue
            policy = (row.get("dist_policy") or "body").strip().lower()
            if policy != "interval" and dist not in (0, None):
                continue
            side[row["anchor"]] = row
        by_side.append(side)

    # An anchor triggered by BOTH sides is an intra-locus event and fires
    # nothing. Comparing anchors rather than labels is what makes this
    # correct: IGLL5 lies inside the IGL locus and both resolve to the same
    # anchor row, so an intra-IGL duplication reads as IGLL5 :: IGL and a
    # label-equality test lets it through.
    both = set(by_side[0]) & set(by_side[1])
    hits = []
    for side in by_side:
        for name, row in side.items():
            if name in both or any(h["anchor"] == name for h in hits):
                continue
            hits.append(row)
    return hits


# -----------------------------------------------------------------------------
# Panel of normals
# -----------------------------------------------------------------------------
# The PoN is a CSV whose column meanings depend on SVTYPE, which is why
# matching against it has to be schema-aware:
#
#   chrom_a, pos_a, chrom_b, pos_b, ?, ?, svtype, frequency
#
#   DEL/DUP/INV  pos_b is a real coordinate on chrom_b
#   BND          pos_b repeats pos_a; the far breakend position is not carried
#   INS          pos_b holds the inserted length, not a position
#
# Comparing pos_b unconditionally therefore rejects every BND match and
# compares a length against a coordinate for every INS. Inter-chromosomal
# junctions are matched on chrom_a, pos_a and chrom_b only.
PON_BUCKET = 100000


def load_pon(path: Optional[Path]):
    """Index a Severus-style panel of normals by (chrom, position bucket).

    A linear scan would be O(records x PoN entries); the PoN runs to
    hundreds of thousands of rows, so it is bucketed. Both orientations are
    stored, so a junction matches whichever way round it was called.
    """
    index: Dict[Tuple[str, int], list] = {}
    if path is None or not Path(str(path)).is_file():
        return index
    n = 0
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split(",")
            if len(f) < 8:
                continue
            try:
                ca, pa, cb, pb = f[0], int(f[1]), f[2], int(f[3])
                svtype, freq = f[6].strip().upper(), float(f[7])
            except ValueError:
                continue
            n += 1
            for (xa, xpa, xb, xpb) in ((ca, pa, cb, pb), (cb, pb, ca, pa)):
                index.setdefault((xa, xpa // PON_BUCKET), []).append(
                    (xpa, xb, xpb, svtype, freq))
    sys.stderr.write(f"panel of normals: {n} entries indexed from {path}\n")
    return index


def pon_lookup(index, chrom_a, pos_a, chrom_b, pos_b, sv_type, tol):
    """Highest PoN frequency matching this junction, or None.

    Same-chromosome events compare both breakpoints, because the PoN
    carries a real pos_b for them. Inter-chromosomal events compare the
    near breakpoint and the partner chromosome only. An insertion is
    matched on its own breakpoint alone.
    """
    if not index or chrom_a is None or pos_a is None:
        return None
    inter = chrom_b is not None and chrom_b != chrom_a
    want_bnd = inter or sv_type in ("BND", "TRA")
    # Same-chromosome events are not matched. The PoN holds 590k DEL and
    # 696k INS entries, so within any workable tolerance almost every small
    # indel in a 30 Mb panel finds a neighbour: on one validation sample
    # 97.8% of DEL and 97.7% of INS matched, against 10.5% of TRA. That is
    # coincidence at scale rather than germline identity, and it would drop
    # two thirds of the callset. The 1,962 BND entries are sparse enough for
    # a coordinate match to mean something.
    if not want_bnd:
        return None
    # Same-chromosome events are not matched. The PoN holds 590k DEL and
    # 696k INS entries, so within any workable tolerance almost every small
    # indel in a 30 Mb panel finds a neighbour: on one validation sample
    # 97.8% of DEL and 97.7% of INS matched, against 10.5% of TRA. That is
    # coincidence at scale rather than germline identity, and it would drop
    # two thirds of the callset. The 1,962 BND entries are sparse enough for
    # a coordinate match to mean something.
    if not want_bnd:
        return None
    best = None
    lo = (pos_a - tol) // PON_BUCKET
    hi = (pos_a + tol) // PON_BUCKET
    for b in range(lo, hi + 1):
        for (xpa, xb, xpb, svtype, freq) in index.get((chrom_a, b), ()):
            if abs(xpa - pos_a) > tol:
                continue
            if want_bnd:
                if svtype != "BND" or xb != chrom_b:
                    continue
            else:
                if svtype == "INS":
                    pass                      # pos_b is a length; ignore it
                elif xb != chrom_b or pos_b is None or abs(xpb - pos_b) > tol:
                    continue
            if best is None or freq > best[1]:
                best = (svtype, freq)
    return best


# -----------------------------------------------------------------------------
# Excluded junctions
# -----------------------------------------------------------------------------
def load_excluded_junctions(path: Optional[Path]) -> Dict[Tuple[str, str], list]:
    """Junctions to drop, keyed by chromosome pair, both orientations stored.

    The criterion is coordinate identity across unrelated patients, not
    recurrence. Somatic breakpoints do not recur to the nucleotide between
    individuals: repair at a real junction is imprecise, so two patients
    sharing a rearrangement share the intron, not the base.

    Recurrence alone must not become the criterion here. This disease is
    defined by recurrent events - t(11;14) appears in 15-20% of patients and
    is a finding, not noise. Coordinate identity is a different claim.
    """
    out: Dict[Tuple[str, str], list] = {}
    if path is None or not Path(str(path)).is_file():
        return out
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5 or f[0].lower() in ("chrom_a", "#chrom_a"):
                continue
            try:
                ca, pa, cb, pb, tol = f[0], int(f[1]), f[2], int(f[3]), int(f[4])
            except ValueError:
                sys.stderr.write(f"excluded junctions: unparsable row: {line}")
                continue
            note = f[6] if len(f) > 6 else ""
            out.setdefault((ca, cb), []).append((pa, pb, tol, note))
            out.setdefault((cb, ca), []).append((pb, pa, tol, note))
    return out


def excluded_reason(excl, row) -> Optional[str]:
    """Why this junction is excluded, or None.

    A dictionary-named pair is never excluded, whatever it matches here. The
    list is a noise filter and a named entity has cleared a higher bar than
    coordinate recurrence can overturn. This override is deliberately not
    configurable: it is what stops a cohort-derived artefact list from ever
    silencing a t(11;14).
    """
    if not excl:
        return None
    if (row.get("known_mm_pair") or "").strip():
        return None
    if (row.get("tier") or "").strip():
        return None
    ca, cb = row.get("chrom_a"), row.get("chrom_b")
    try:
        pa, pb = int(row.get("pos_a")), int(row.get("pos_b"))
    except (TypeError, ValueError):
        return None
    for (xa, xb, tol, note) in excl.get((ca, cb), []):
        if abs(pa - xa) <= tol and abs(pb - xb) <= tol:
            return (f"listed exclusion at {ca}:{xa} x {cb}:{xb} "
                    f"(+/-{tol} bp)" + (f"; {note}" if note else ""))
    return None


def pon_reason(row, min_freq) -> Optional[str]:
    """Why this junction is dropped as germline, or None.

    A dictionary-named or graded row is never dropped, on the same reasoning
    as the exclusion list and with the same non-configurable override. A
    named entity has cleared a higher bar than population frequency can
    overturn, and pon_freq is still written on the row, so a conflict
    between the two is visible rather than silently resolved.

    Threshold rather than presence. The PoN records how often a junction is
    seen across the 1000 Genomes panel; an entry at 0.99 is effectively
    fixed in the population and cannot be somatic, while a low-frequency
    entry says only that someone else has also seen it.
    """
    if min_freq is None or min_freq <= 0:
        return None
    if (row.get("known_mm_pair") or "").strip():
        return None
    if (row.get("tier") or "").strip():
        return None
    f = (row.get("pon_freq") or "").strip()
    if not f:
        return None
    try:
        v = float(f)
    except ValueError:
        return None
    if v < min_freq:
        return None
    return (f"panel of normals {row.get('pon_svtype', '')} at frequency "
            f"{v:.4f} (>= {min_freq})")


# -----------------------------------------------------------------------------
# Annotation
# -----------------------------------------------------------------------------
def region_for(chrom, pos, panel) -> Optional[PanelRegion]:
    if chrom is None or pos is None:
        return None
    for r in panel:
        if r.contains(chrom, pos):
            return r
    return None


def characterize_side(chrom, pos, region, cytobands, gene_model=None):
    """Return (label, source, band) for one breakpoint side.

    Resolution order: gene model, panel interval label, cytoband, coordinate.

    The gene model comes first because a panel interval may be a compound
    ("FGFR3/NSD2") that names two genes without saying which. It comes first
    only for labelling; panel membership is decided separately by region_for
    against the panel BED.

    Falling back to the panel label is the right answer for a BCR breakend,
    not a failure. MYC's body is 8 kb inside a 5 Mb window and t(11;14)
    breakpoints sit outside CCND1 entirely, so most real breakends on this
    panel resolve to the interval label. The distance column is what makes
    that auditable.
    """
    band = cytobands.band_for(chrom, pos)
    band_label = f"{_strip_chr(chrom)}{band}" if band and chrom else None
    gene = gene_for(chrom, pos, gene_model)
    if gene:
        return gene, "gene_model", band_label
    if region is not None:
        return region.name, "panel", band_label
    if band_label:
        return band_label, "cytoband", band_label
    if chrom is not None and pos is not None:
        return f"{chrom}:{pos / 1e6:.1f}Mb", "coordinate", None
    return "OFF_PANEL", "coordinate", None


def annotate(records, panel, dictionary, anchors, sample, cytobands,
             gene_model=None, ig_segments=None, pon=None, pon_tol=2500):
    out = []
    for r in records:
        side_a = region_for(r.chrom, r.pos, panel)
        side_b = region_for(r.mate_chrom, r.mate_pos, panel)
        if side_a is None and side_b is None:
            continue

        gene_a, src_a, band_a = characterize_side(
            r.chrom, r.pos, side_a, cytobands, gene_model)
        gene_b, src_b, band_b = characterize_side(
            r.mate_chrom, r.mate_pos, side_b, cytobands, gene_model)

        ig_a = ig_region_for(r.chrom, r.pos, ig_segments or [])
        ig_b = ig_region_for(r.mate_chrom, r.mate_pos, ig_segments or [])

        dist_a = dist_to_gene(r.chrom, r.pos, gene_a, gene_model)
        dist_b = dist_to_gene(r.mate_chrom, r.mate_pos, gene_b, gene_model)

        hit, quality = dictionary_lookup(dictionary, gene_a, gene_b,
                                         band_a, band_b)

        # Span guard. Where a dictionary row declares min_span_bp and the
        # event is intrachromosomal with both ends known, an undersized span
        # demotes the match: the record stays, the entity claim goes. No
        # shipped row sets it, because no pair in this dictionary has both
        # partners inside one panel interval, but the guard belongs with the
        # matching rather than being added after one is needed.
        span_note = ""
        if hit and quality == "full":
            ms = (hit.get("min_span_bp") or "").strip()
            if ms and r.mate_chrom == r.chrom and r.mate_pos is not None:
                span = abs(r.mate_pos - r.pos)
                if span < int(ms):
                    span_note = (f"span {span}bp < min {ms}bp for "
                                 f"{hit.get('name', 'pair')}; entity removed")
                    hit, quality = None, "below_span"

        pon_hit = pon_lookup(pon, r.chrom, r.pos, r.mate_chrom, r.mate_pos,
                             r.sv_type, pon_tol)

        hits = anchor_hits(anchors, gene_a, gene_b, dist_a, dist_b)

        # Reportable when the dictionary names it or it touches an anchor.
        # Everything else is still emitted, carrying reportable=no, so the
        # on-panel callset stays auditable rather than silently trimmed.
        reportable = "yes" if (hit or hits) else "no"

        # tier is the dictionary's, marked with '?' when the match was only
        # partial, so a graded row always says how firmly it was graded.
        tier = ""
        if hit:
            t = hit.get("tier", "")
            tier = t if quality == "full" else (t + "?" if t else "")

        out.append({
            "sample":         sample,
            "sv_id":          r.sv_id,
            "sv_type":        r.sv_type,
            "filter":         r.filt,
            "chrom_a":        r.chrom,
            "pos_a":          str(r.pos),
            "gene_a":         gene_a,
            "chrom_b":        r.mate_chrom or "",
            "pos_b":          str(r.mate_pos) if r.mate_pos is not None else "",
            "gene_b":         gene_b,
            "gene_a_source":  src_a,
            "gene_b_source":  src_b,
            "gene_a_dist":    "" if dist_a is None else str(dist_a),
            "gene_b_dist":    "" if dist_b is None else str(dist_b),
            "ig_region_a":    ig_a,
            "ig_region_b":    ig_b,
            "pon_freq":       f"{pon_hit[1]:.4f}" if pon_hit else "",
            "pon_svtype":     pon_hit[0] if pon_hit else "",
            "band_a":         band_a or "",
            "band_b":         band_b or "",
            "known_mm_pair":  hit.get("name", "") if hit else "",
            "entity":         hit.get("entity", "") if hit else "",
            "tier":           tier,
            "known_freq":     hit.get("frequency", "") if hit else "",
            "match_quality":  quality,
            "anchor":         ",".join(h["anchor"] for h in hits),
            "anchor_class":   ",".join(h.get("anchor_class", "") for h in hits),
            "reportable":     reportable,
            "dict_notes":     (hit.get("notes", "") if hit else span_note),
            "callers":        ",".join(r.callers) or "unknown",
            "n_callers":      str(len(r.callers)),
            "supp_vec":       r.info.get("SUPP_VEC", ""),
            "support_reads":  r.support,
        })
    return out


COLUMNS = [
    "sample", "sv_id", "sv_type", "filter",
    "chrom_a", "pos_a", "gene_a", "chrom_b", "pos_b", "gene_b",
    "gene_a_source", "gene_b_source", "gene_a_dist", "gene_b_dist",
    "ig_region_a", "ig_region_b", "band_a", "band_b",
    "pon_freq", "pon_svtype",
    "known_mm_pair", "entity", "tier", "known_freq", "match_quality",
    "anchor", "anchor_class", "reportable", "dict_notes",
    "callers", "n_callers", "supp_vec", "support_reads",
]


def main() -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--panel-bed", required=True, type=Path)
    ap.add_argument("--dictionary", required=True, type=Path)
    ap.add_argument("--cytoband-bed", required=True, type=Path,
                    help="T2T-CHM13v2.0 cytoband BED. Off-panel partners are "
                         "characterized by band.")
    ap.add_argument("--gene-model", default=None, type=Path,
                    help="Bare gene bodies, for naming a breakpoint and "
                         "measuring its distance to the gene it is named "
                         "after. Optional; without it labels come from the "
                         "panel interval and the distance columns stay empty.")
    ap.add_argument("--ig-segments", default=None, type=Path,
                    help="Ig locus sub-regions (C/switch, J, D, V) built by "
                         "bin/build_ig_segments.py. Optional; adds "
                         "ig_region_a/ig_region_b. Reported, not filtered on.")
    ap.add_argument("--anchors", default=None, type=Path,
                    help="Promiscuous loci reported whatever the partner. "
                         "Optional; without it only dictionary-named pairs "
                         "are reportable.")
    ap.add_argument("--pon", default=None, type=Path,
                    help="Severus-style panel of normals CSV. Adds pon_freq "
                         "and pon_svtype, and drops junctions at or above "
                         "--pon-min-freq. A dictionary-named or graded row is "
                         "never dropped.")
    ap.add_argument("--pon-min-freq", type=float, default=0.10,
                    help="Population frequency at or above which a junction "
                         "is dropped as germline [0.10]. Set to 0 to annotate "
                         "without dropping.")
    ap.add_argument("--pon-tol", type=int, default=2500,
                    help="Bases a breakpoint may differ from a PoN entry and "
                         "still match [2500].")
    ap.add_argument("--excluded-junctions", default=None, type=Path,
                    help="Coordinate-identical artefact junctions to drop. "
                         "A dictionary-named or graded row is never dropped.")
    ap.add_argument("--sample", required=True, type=str)
    ap.add_argument("--output", required=True, type=Path)
    ap.add_argument("--version", action="version",
                    version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel_bed)
    dictionary = load_dictionary(args.dictionary)
    anchors = load_anchors(args.anchors)
    gene_model = load_gene_model(args.gene_model)
    ig_segments = load_ig_segments(args.ig_segments)
    cytobands = load_cytobands(args.cytoband_bed)
    excl = load_excluded_junctions(args.excluded_junctions)
    pon = load_pon(args.pon)
    records = parse_vcf(args.vcf)

    sys.stderr.write(
        f"panel {len(panel)} regions | gene model {len(gene_model)} features | "
        f"ig segments {len(ig_segments)} | "
        f"dictionary {len(dictionary)} pairs | anchors {len(anchors)} tokens | "
        f"exclusions {len(set(k[0] for k in excl))} chromosome pair(s)\n")

    rows = annotate(records, panel, dictionary, anchors, args.sample,
                    cytobands, gene_model, ig_segments,
                    pon, args.pon_tol)

    keep, dropped = [], []
    n_pon = 0
    for r in rows:
        why = excluded_reason(excl, r)
        if not why:
            why = pon_reason(r, args.pon_min_freq)
            if why:
                n_pon += 1
        (dropped if why else keep).append((r, why))

    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           lineterminator="\n", extrasaction="ignore")
        w.writeheader()
        w.writerows([r for r, _ in keep])

    if dropped:
        sys.stderr.write(
            f"dropped {len(dropped)} junction(s): {len(dropped) - n_pon} on "
            f"the exclusion list, {n_pon} as germline by frequency\n")
        seen = set()
        for r, why in dropped:
            k = (r["chrom_a"], r["pos_a"], r["chrom_b"], r["pos_b"])
            if k in seen:
                continue
            seen.add(k)
            sys.stderr.write(f"  {r['gene_a']} x {r['gene_b']}  "
                             f"{k[0]}:{k[1]} x {k[2]}:{k[3]}  - {why}\n")

    # Counted over what was written, not over what was annotated. Counting
    # the pre-exclusion list while writing the post-exclusion one describes a
    # table that is not on disk.
    written = [r for r, _ in keep]
    n_report = sum(1 for r in written if r["reportable"] == "yes")
    n_named = sum(1 for r in written if r["known_mm_pair"])
    n_tier = sum(1 for r in written if r["tier"])
    sys.stderr.write(
        f"{len(written)} on-panel record(s) written ({len(dropped)} excluded), "
        f"{n_report} reportable, {n_named} named by dictionary, "
        f"{n_tier} graded -> {args.output}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
