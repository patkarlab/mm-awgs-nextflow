#!/usr/bin/env python3
"""
annotate_mm_translocations.py
=============================

Annotates a SURVIVOR-merged VCF of structural variants against the MM-specific
translocation panel by intersecting BND breakpoints with the panel BED.

For each record where at least one breakpoint falls inside a panel region,
emits a row of:
  - sample id
  - SV id / type / filter
  - both breakpoints (chrom, pos, gene region)
  - whether the pair matches a known MM partner pair (from dictionary)
  - supporting callers inferred from SURVIVOR's SUPP_VEC
  - support read counts from RE / SR

The dictionary file is the only source of biological priors. The script
itself never hardcodes any expected breakpoints, sample-specific findings,
or known karyotypes.

Usage:
  annotate_mm_translocations.py \\
      --vcf merged.vcf.gz \\
      --panel-bed aWGS_MMfocused_v6_t2t_chr.bed \\
      --dictionary mm_translocation_dictionary.tsv \\
      --sample 11F20262905 \\
      --output 11F20262905.mm_annotated.tsv
"""

from __future__ import annotations

import argparse
import gzip
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


__version__ = "0.1.0"


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


def load_dictionary(dict_path: Path) -> Dict[Tuple[str, str], Dict[str, str]]:
    """
    Load the MM translocation dictionary, keyed by unordered (a, b) pair
    of uppercase symbols. Missing dictionary file is non-fatal.
    """
    out: Dict[Tuple[str, str], Dict[str, str]] = {}
    if not dict_path.exists():
        return out
    with open(dict_path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            a = (row.get("partner_a") or "").strip().upper()
            b = (row.get("partner_b") or "").strip().upper()
            if not a or not b:
                continue
            key = tuple(sorted([a, b]))
            out[key] = row
    return out


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
            out.append(SvRecord(
                chrom=chrom, pos=pos, sv_id=sv_id, sv_type=sv_type,
                mate_chrom=mate_chrom, mate_pos=mate_pos,
                filt=filt, info=info, callers=infer_callers(info),
            ))
    return out


def region_for(chrom: Optional[str], pos: Optional[int], panel: List[PanelRegion]) -> Optional[PanelRegion]:
    if chrom is None or pos is None:
        return None
    for r in panel:
        if r.contains(chrom, pos):
            return r
    return None


def annotate(records, panel, dictionary, sample):
    out = []
    for r in records:
        side_a = region_for(r.chrom, r.pos, panel)
        side_b = region_for(r.mate_chrom, r.mate_pos, panel)
        if side_a is None and side_b is None:
            continue

        gene_a = side_a.name if side_a else "OFF_PANEL"
        gene_b = side_b.name if side_b else "OFF_PANEL"

        known = ""
        freq = ""
        if side_a and side_b:
            key = tuple(sorted([gene_a.upper(), gene_b.upper()]))
            hit = dictionary.get(key)
            if hit:
                known = hit.get("name", "yes")
                freq = hit.get("frequency", "")

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
            "known_mm_pair":  known,
            "known_freq":     freq,
            "callers":        ",".join(r.callers) or "unknown",
            "n_callers":      str(len(r.callers)),
            "supp_vec":       r.info.get("SUPP_VEC", ""),
            "support_reads":  r.info.get("RE", r.info.get("SR", "")),
        })
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf",         required=True, type=Path)
    ap.add_argument("--panel-bed",   required=True, type=Path)
    ap.add_argument("--dictionary",  required=True, type=Path)
    ap.add_argument("--sample",      required=True, type=str)
    ap.add_argument("--output",      required=True, type=Path)
    ap.add_argument("--version",     action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel_bed)
    dictionary = load_dictionary(args.dictionary)
    records = parse_vcf(args.vcf)
    rows = annotate(records, panel, dictionary, args.sample)

    columns = [
        "sample", "sv_id", "sv_type", "filter",
        "chrom_a", "pos_a", "gene_a",
        "chrom_b", "pos_b", "gene_b",
        "known_mm_pair", "known_freq",
        "callers", "n_callers", "supp_vec", "support_reads",
    ]
    with open(args.output, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row.get(c, "") for c in columns) + "\n")

    print(f"Annotated {len(rows)} on-panel SV records -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
