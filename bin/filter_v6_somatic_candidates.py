#!/usr/bin/env python3
"""
filter_v6_somatic_candidates.py

Reduce a Clair3 + VEP somatic-candidate TSV to a clinically reviewable report
for the adaptive-WGS multiple myeloma project (v6 panel).

Standard-library only (csv) - no pandas/numpy - so it runs in any python3,
including the awgs_sv env, and as a Nextflow bin/ script.

Hard filters (remove rows):
  1. PANEL-GENE membership : keep variants whose gene symbol is in the v6 panel.
     Keyed on gene SYMBOL, not coordinates, because candidate TSVs are hg38
     while the v6 BED is T2T - symbols are build-independent.
  2. CONSEQUENCE class      : keep protein-altering consequences (the set
     selected in the Excel screenshot), matched on the most-severe (first) token.

Soft filter (labels only, never removes):
  COMMON_POLYMORPHISM : pop_af_max > --max-pop-af (default 0.01), where present.

Column hygiene applied to every output:
  - '.'/'' -> '-1' sentinel across all fields.
  - exon 'rank/total' split into exon_rank/exon_total integers (drops the
    Excel date-corrupting 'exon' field).
  - tumor_af_pct added beside tumor_af.
  - domains dropped.
  - REF_COUNT/ALT_COUNT/DP passed through if present (populated upstream by the
    VEP step or augment_alt_counts.py), else emitted as -1.

No alt-read-support filter is applied; REF_COUNT/ALT_COUNT/DP are reported.

Outputs per input <stem>:
  <stem>.v6_filtered.tsv   in-panel + reportable rows, with Filter column
  <stem>.v6_clinical.tsv   Filter == PASS only, sorted by gene then position
"""

import argparse
import csv
import os
import sys

MISSING = "-1"

PANEL_GENE_TO_LABEL = {
    "CDKN2C": "CDKN2C",
    "NRAS": "NRAS",
    "TENT5C": "TENT5C", "FAM46C": "TENT5C",
    "FCRL4": "FCRL4_FCRL5", "FCRL5": "FCRL4_FCRL5",
    "CXCR4": "CXCR4",
    "ATR": "ATR",
    "FGFR3": "FGFR3_NSD2",
    "NSD2": "FGFR3_NSD2", "WHSC1": "FGFR3_NSD2", "MMSET": "FGFR3_NSD2",
    "TNFAIP8": "TNFAIP8",
    "EGR1": "EGR1",
    "IRF4": "IRF4",
    "TXNDC5": "TXNDC5",
    "H1-4": "H1-4", "HIST1H1E": "H1-4", "H1F4": "H1-4",
    "LTB": "LTB",
    "CCND3": "CCND3",
    "PRDM1": "PRDM1",
    "BRAF": "BRAF",
    "MYC": "MYC",
    "MAFA": "MAFA",
    "PAX5": "PAX5",
    "CCND1": "CCND1",
    "ATM": "ATM",
    "CCND2": "CCND2",
    "KRAS": "KRAS",
    "LRRK2": "LRRK2",
    "RB1": "RB1",
    "DIS3": "DIS3",
    "MAX": "MAX",
    "TRAF3": "TRAF3",
    "CYLD": "CYLD",
    "WWOX": "WWOX_MAF", "MAF": "WWOX_MAF",
    "TP53": "TP53",
    "MAP3K14": "MAP3K14", "NIK": "MAP3K14",
    "BCL2": "BCL2",
    "MAFB": "MAFB",
    "XBP1": "XBP1",
}

IG_PREFIX_TO_LABEL = {"IGH": "IGH_locus", "IGK": "IGK_locus", "IGL": "IGL_locus"}

REPORTABLE_CONSEQUENCES = {
    "frameshift_variant",
    "inframe_deletion",
    "inframe_insertion",
    "missense_variant",
    "start_lost",
    "stop_gained",
    "stop_lost",
    "stop_retained_variant",
}

# Final column order; any unexpected columns are appended after these.
PREFERRED_COLS = [
    "chrom", "pos", "ref", "alt", "qual", "variant_type",
    "gene", "panel_label", "transcript", "biotype", "canonical",
    "consequence", "csq_primary", "impact", "exon_rank", "exon_total",
    "rs_id", "pop_af_max", "pop_af_max_source", "clinvar_sig",
    "tumor_af", "tumor_af_pct", "REF_COUNT", "ALT_COUNT", "DP", "Filter",
]

MISSING_TOKENS = {"", "-1", ".", "nan", "NA"}


def is_missing(value):
    return value is None or str(value).strip() in MISSING_TOKENS


def safe_float(value):
    """Return float(value) or None if missing/non-numeric."""
    if is_missing(value):
        return None
    try:
        return float(value)
    except (ValueError, TypeError):
        return None


def clean_cell(value):
    if value is None:
        return MISSING
    s = str(value).strip()
    return MISSING if s in ("", ".") else s


def resolve_panel_label(gene_symbol, include_ig):
    g = str(gene_symbol).strip()
    if g in PANEL_GENE_TO_LABEL:
        return PANEL_GENE_TO_LABEL[g]
    if include_ig:
        for prefix, label in IG_PREFIX_TO_LABEL.items():
            if g.upper().startswith(prefix):
                return label
    return None


def split_exon(value):
    s = str(value).strip()
    if s in ("", "-1", "."):
        return MISSING, MISSING
    if "/" in s:
        rank, total = s.split("/", 1)
        return (rank.strip() or MISSING), (total.strip() or MISSING)
    return s, MISSING


def af_to_pct(value):
    f = safe_float(value)
    if f is None:
        return MISSING
    return f"{f * 100:.2f}"


def filter_one_file(input_path, outdir, max_pop_af, include_ig):
    stem = os.path.basename(input_path)
    for suffix in (".tsv", ".txt"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]

    with open(input_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        in_cols = list(reader.fieldnames or [])
        rows = [dict(r) for r in reader]
    n_input = len(rows)

    kept = []
    for row in rows:
        # 1. Normalise missing values.
        for key in list(row.keys()):
            row[key] = clean_cell(row[key])

        # 2. Ensure count columns exist.
        for key in ("REF_COUNT", "ALT_COUNT", "DP"):
            row.setdefault(key, MISSING)
            if is_missing(row[key]):
                row[key] = MISSING

        # 2a. Exon split (and drop original 'exon').
        if "exon" in row:
            row["exon_rank"], row["exon_total"] = split_exon(row.pop("exon"))

        # 2b. AF percentage.
        if "tumor_af" in row:
            row["tumor_af_pct"] = af_to_pct(row["tumor_af"])

        # 2c. Drop domains.
        row.pop("domains", None)

        # 3. Panel membership + primary consequence.
        label = resolve_panel_label(row.get("gene", ""), include_ig)
        primary = str(row.get("consequence", "")).strip().split("&")[0]

        # 4. Hard filters.
        if label is None or primary not in REPORTABLE_CONSEQUENCES:
            continue
        row["panel_label"] = label
        row["csq_primary"] = primary

        # 5. Soft polymorphism flag.
        pop_af = safe_float(row.get("pop_af_max"))
        row["Filter"] = "COMMON_POLYMORPHISM" if (pop_af is not None and pop_af > max_pop_af) else "PASS"

        kept.append(row)

    # Sort by panel gene, chrom, numeric position.
    def sort_key(r):
        return (r.get("panel_label", ""), r.get("chrom", ""), safe_float(r.get("pos")) or 0.0)
    kept.sort(key=sort_key)

    # Column order.
    all_cols = set()
    for r in kept:
        all_cols.update(r.keys())
    ordered = [c for c in PREFERRED_COLS if c in all_cols]
    ordered += [c for c in all_cols if c not in ordered]

    os.makedirs(outdir, exist_ok=True)
    filtered_path = os.path.join(outdir, f"{stem}.v6_filtered.tsv")
    clinical_path = os.path.join(outdir, f"{stem}.v6_clinical.tsv")

    def write_tsv(path, records):
        with open(path, "w", newline="") as out:
            w = csv.DictWriter(out, fieldnames=ordered, delimiter="\t", lineterminator="\n",
                               extrasaction="ignore", restval=MISSING)
            w.writeheader()
            for rec in records:
                w.writerow(rec)

    write_tsv(filtered_path, kept)
    clinical = [r for r in kept if r["Filter"] == "PASS"]
    write_tsv(clinical_path, clinical)

    n_pass = len(clinical)
    n_common = sum(1 for r in kept if r["Filter"] == "COMMON_POLYMORPHISM")
    genes_hit = sorted({r["panel_label"] for r in clinical})
    sys.stderr.write(
        f"[{stem}] input={n_input}  in-panel+reportable={len(kept)}  "
        f"PASS={n_pass}  COMMON_POLYMORPHISM={n_common}\n"
    )
    sys.stderr.write(f"[{stem}] PASS genes: {', '.join(genes_hit) if genes_hit else '(none)'}\n")
    return {"sample": stem, "input": n_input, "kept": len(kept),
            "pass": n_pass, "common_polymorphism": n_common}


def parse_args():
    ap = argparse.ArgumentParser(
        description="Filter Clair3/VEP somatic candidates to the v6 MM panel.")
    ap.add_argument("-i", "--input", required=True, nargs="+",
                    help="One or more *_somatic_candidates.tsv files.")
    ap.add_argument("-o", "--outdir", default="v6_filtered",
                    help="Output directory (default: ./v6_filtered).")
    ap.add_argument("--max-pop-af", type=float, default=0.01,
                    help="pop_af_max above this is flagged COMMON_POLYMORPHISM "
                         "(default: 0.01).")
    ap.add_argument("--include-ig", action="store_true",
                    help="Include IGH/IGK/IGL-locus genes (off by default: "
                         "IG-locus SNVs are dominated by somatic hypermutation).")
    return ap.parse_args()


def main():
    args = parse_args()
    summary_rows = []
    for path in args.input:
        if not os.path.isfile(path):
            sys.stderr.write(f"ERROR: input not found: {path}\n")
            sys.exit(1)
        summary_rows.append(filter_one_file(path, args.outdir, args.max_pop_af, args.include_ig))

    os.makedirs(args.outdir, exist_ok=True)
    summary_path = os.path.join(args.outdir, "v6_filter_summary.tsv")
    with open(summary_path, "w", newline="") as out:
        cols = ["sample", "input", "kept", "pass", "common_polymorphism"]
        w = csv.DictWriter(out, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in summary_rows:
            w.writerow(r)
    sys.stderr.write(f"\nWrote cohort summary: {summary_path}\n")


if __name__ == "__main__":
    main()
