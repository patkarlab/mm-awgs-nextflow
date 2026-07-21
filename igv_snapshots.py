#!/usr/bin/env python3
"""
igv_snapshots.py

Generate self-contained igv-reports HTML snapshots for a single sample, for
either of two evidence classes in the MM adaptive-WGS pipeline:

  somatic         one IGV view per v6 clinical somatic SNV/indel, rendered
                  against the hg38 BAM. Input is the v6 clinical TSV
                  (Filter == PASS, in-panel) produced by
                  bin/filter_v6_somatic_candidates.py.

  translocations  one IGV view per breakpoint, rendered against the T2T BAM.
                  Each translocation row carries two breakpoints (A and B);
                  this script explodes each row into two single-locus sites
                  sharing one EVENT id, so both ends of a BND are inspectable
                  and stay linked in the variant table. Input is the
                  <sample>.translocations.tsv produced by merge_translocations.py.

Design notes
------------
- Standard library only (csv, argparse, subprocess, os, sys, tempfile). Runs in
  the awgs_sv env, which has no pandas. igv-reports' `create_report` must be on
  PATH (it is, in awgs_sv: igv_reports 0.1.0).
- create_report is driven through a tab-delimited *sites* file (not a VCF), so
  the variant table shows exactly the columns we write into that sites file.
  Column-to-coordinate mapping is passed via --sequence/--begin/--end.
- This script never hardcodes any variant, gene-pair, FISH finding, or expected
  call. It reads whatever rows the upstream filters produced and renders them.
  The only gene-aware content is column *names* read from the input header.
- Output is one HTML per sample per mode. Empty input (zero rows) is a normal,
  expected outcome (a sample may have no clinical somatic variants, or no
  annotated translocations); the script writes a small placeholder HTML and
  exits 0 so it composes cleanly into a per-sample report.

Exit codes
----------
0  success, including the legitimate "no sites" case
2  usage / input error (missing BAM, missing input TSV, create_report failure)
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def eprint(*args):
    """Print to stderr (stdout is reserved for any machine-readable summary)."""
    print(*args, file=sys.stderr)


def fail(msg, code=2):
    eprint("ERROR: " + msg)
    sys.exit(code)


def read_tsv(path):
    """Read a tab-delimited file with a header row into (header, rows).

    rows is a list of dicts keyed by the header names. Returns ([], []) if the
    file is empty or header-only.
    """
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return [], []
        rows = []
        for raw in reader:
            if not raw:
                continue
            # Pad short rows so dict zip does not silently drop trailing cols.
            if len(raw) < len(header):
                raw = raw + [""] * (len(header) - len(raw))
            rows.append(dict(zip(header, raw)))
        return header, rows


def write_placeholder_html(out_html, title, message):
    """Write a minimal standalone HTML for the legitimate 'no sites' case.

    Kept deliberately dependency-free so it renders anywhere and slots into the
    dashboard the same way a real report would.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    html = (
        "<!DOCTYPE html><html><head><meta charset='utf-8'>"
        "<title>{title}</title></head><body>"
        "<h2>{title}</h2><p>{message}</p>"
        "</body></html>"
    ).format(title=title, message=message)
    with open(out_html, "w") as fh:
        fh.write(html)


def index_exists_and_fresh(bam):
    """Return True iff a .bai exists and is at least as new as the BAM.

    The pipeline has previously hit stale .bai files (BAM rewritten after
    indexing) causing mosdepth/igv-reports BGZF decode errors. We treat an
    older index as absent so the caller can re-index.
    """
    bai = bam + ".bai"
    alt = bam[:-4] + ".bai" if bam.endswith(".bam") else None
    candidates = [p for p in (bai, alt) if p]
    for p in candidates:
        if os.path.exists(p) and os.path.getmtime(p) >= os.path.getmtime(bam):
            return True
    return False


def ensure_bam_index(bam, samtools="samtools"):
    """Index the BAM if no fresh index is present. Safe to call repeatedly."""
    if index_exists_and_fresh(bam):
        return
    eprint("[index] (re)indexing {0}".format(bam))
    try:
        subprocess.run([samtools, "index", bam], check=True)
    except subprocess.CalledProcessError as exc:
        fail("samtools index failed on {0}: {1}".format(bam, exc))
    except FileNotFoundError:
        fail("samtools not found on PATH; cannot index {0}".format(bam))


# ---------------------------------------------------------------------------
# Mode: somatic
# ---------------------------------------------------------------------------

# Columns we want in the variant table, in display order, IF present in the
# input header. We never require a column; whatever exists is shown.
SOMATIC_TABLE_COLUMNS = [
    "gene", "panel_label", "consequence", "impact",
    "tumor_af_pct", "ALT_COUNT", "DP",
    "clinvar_sig", "rs_id", "Filter",
]

# Required coordinate columns for the somatic TSV.
SOMATIC_CHROM_COL = "chrom"
SOMATIC_POS_COL = "pos"
SOMATIC_REF_COL = "ref"
SOMATIC_ALT_COL = "alt"


def build_somatic_sites(rows, header, sites_path):
    """Write the igv-reports sites file for somatic variants.

    The sites file is the v6 clinical rows, re-emitted with a guaranteed
    coordinate layout: chrom, start, end, then the descriptive columns that
    exist. start/end are 1-based identical (point variants); igv-reports'
    --flanking widens the view. For indels we still anchor on POS; the BAM
    pileup shows the event.
    """
    table_cols = [c for c in SOMATIC_TABLE_COLUMNS if c in header]
    # Front-load ref/alt right after coords if available — useful in the table.
    lead_cols = [c for c in (SOMATIC_REF_COL, SOMATIC_ALT_COL) if c in header]
    out_header = ["chrom", "start", "end"] + lead_cols + table_cols

    n_written = 0
    with open(sites_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(out_header)
        for r in rows:
            chrom = r.get(SOMATIC_CHROM_COL, "").strip()
            pos = r.get(SOMATIC_POS_COL, "").strip()
            if not chrom or not pos:
                continue
            try:
                pos_i = int(pos)
            except ValueError:
                continue
            line = [chrom, str(pos_i), str(pos_i)]
            line += [r.get(c, "") for c in lead_cols]
            line += [r.get(c, "") for c in table_cols]
            writer.writerow(line)
            n_written += 1
    return out_header, n_written


# ---------------------------------------------------------------------------
# Mode: translocations
# ---------------------------------------------------------------------------

# Coordinate columns in the translocations TSV (from merge_translocations.py).
TRA_CHROM_A = "chrom_a"
TRA_POS_A = "pos_a"
TRA_GENE_A = "gene_a"
TRA_CHROM_B = "chrom_b"
TRA_POS_B = "pos_b"
TRA_GENE_B = "gene_b"

# Descriptive columns to carry into the table for each breakpoint, if present.
TRA_TABLE_COLUMNS = [
    "event", "side", "partner_locus",
    "sv_id", "sv_type", "filter",
    "known_mm_pair", "known_freq",
    "callers", "n_callers", "support_reads",
    "support_sniffles", "support_cutesv", "support_severus",
]


def build_translocation_sites(rows, header, sites_path):
    """Write the igv-reports sites file, two breakpoints per translocation row.

    Each input row becomes two output rows (side A, side B). Both share an
    'event' id (the row's sv_id, or a synthesised index) and each records its
    partner locus, so the two ends remain associated in the table even though
    igv-reports renders each as its own single-locus view.
    """
    have_a = TRA_CHROM_A in header and TRA_POS_A in header
    have_b = TRA_CHROM_B in header and TRA_POS_B in header
    if not have_a:
        fail("translocations TSV lacks {0}/{1} columns"
             .format(TRA_CHROM_A, TRA_POS_A))

    # Descriptive columns actually present (minus the synthesised ones we add).
    synth = {"event", "side", "partner_locus"}
    passthrough = [c for c in TRA_TABLE_COLUMNS
                   if c not in synth and c in header]
    out_header = (["chrom", "start", "end", "event", "side", "partner_locus"]
                  + passthrough)

    def emit(writer, chrom, pos, event, side, partner, src_row):
        try:
            pos_i = int(str(pos).strip())
        except (ValueError, AttributeError):
            return 0
        line = [chrom, str(pos_i), str(pos_i), event, side, partner]
        line += [src_row.get(c, "") for c in passthrough]
        writer.writerow(line)
        return 1

    n_written = 0
    with open(sites_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(out_header)
        for idx, r in enumerate(rows):
            event = (r.get("sv_id") or "event_{0}".format(idx)).strip()
            chrom_a = r.get(TRA_CHROM_A, "").strip()
            pos_a = r.get(TRA_POS_A, "").strip()
            gene_a = r.get(TRA_GENE_A, "").strip()
            chrom_b = r.get(TRA_CHROM_B, "").strip()
            pos_b = r.get(TRA_POS_B, "").strip()
            gene_b = r.get(TRA_GENE_B, "").strip()

            partner_b = "{0}:{1}".format(chrom_b or "?", gene_b or "?")
            partner_a = "{0}:{1}".format(chrom_a or "?", gene_a or "?")

            if chrom_a and pos_a:
                n_written += emit(writer, chrom_a, pos_a, event,
                                  "A:{0}".format(gene_a or "?"),
                                  partner_b, r)
            if have_b and chrom_b and pos_b:
                n_written += emit(writer, chrom_b, pos_b, event,
                                  "B:{0}".format(gene_b or "?"),
                                  partner_a, r)
    return out_header, n_written


# ---------------------------------------------------------------------------
# create_report invocation
# ---------------------------------------------------------------------------

def run_create_report(sites_path, fasta, bam, out_html, flanking,
                      title, create_report="create_report",
                      exclude_flags="1536"):
    """Assemble and run the create_report command for a tab-delimited sites file.

    --sequence/--begin/--end use 1-based column NUMBERS (igv-reports counts the
    sites-file columns from 1). Our sites files always put chrom=1, start=2,
    end=3, so those are fixed.

    exclude_flags default 1536 = 0x600 = unmapped(4)? no: 1024 (PCR/optical
    duplicate) + 512 (fails QC). We deliberately do NOT exclude secondary or
    supplementary alignments, because split/supplementary reads are the
    evidence for a translocation breakpoint and must remain visible.
    """
    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    cmd = [
        create_report,
        sites_path,
        "--fasta", fasta,
        "--sequence", "1",
        "--begin", "2",
        "--end", "3",
        "--flanking", str(flanking),
        "--tracks", bam,
        "--standalone",
        "--exclude-flags", str(exclude_flags),
        "--title", title,
        "--output", out_html,
    ]
    eprint("[create_report] " + " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        fail("create_report failed (exit {0}) for {1}"
             .format(exc.returncode, out_html))
    except FileNotFoundError:
        fail("create_report not found on PATH (expected in awgs_sv)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Render igv-reports HTML snapshots for one sample.")
    p.add_argument("--mode", required=True,
                   choices=["somatic", "translocations"],
                   help="Which evidence class to render.")
    p.add_argument("--sample", required=True,
                   help="Sample id (sequencing id only; used in titles/paths).")
    p.add_argument("--sites-tsv", required=True,
                   help="Input TSV: v6 clinical (somatic) or "
                        "translocations TSV (translocations).")
    p.add_argument("--bam", required=True,
                   help="Alignment to display: hg38 BAM for somatic, "
                        "T2T BAM for translocations.")
    p.add_argument("--fasta", required=True,
                   help="Reference FASTA matching --bam "
                        "(hg38 for somatic, T2T chr-named for translocations).")
    p.add_argument("--out-html", required=True,
                   help="Output HTML path.")
    p.add_argument("--flanking", type=int, default=5000,
                   help="bp of context either side of each site (default 5000).")
    p.add_argument("--samtools", default="samtools",
                   help="samtools executable (default: samtools on PATH).")
    p.add_argument("--create-report", default="create_report",
                   help="create_report executable (default: on PATH).")
    args = p.parse_args(argv)

    # Validate inputs early with clear messages.
    if not os.path.exists(args.sites_tsv):
        fail("input TSV not found: {0}".format(args.sites_tsv))
    if not os.path.exists(args.bam):
        fail("BAM not found: {0}".format(args.bam))
    if not os.path.exists(args.fasta):
        fail("FASTA not found: {0}".format(args.fasta))
    fai = args.fasta + ".fai"
    if not os.path.exists(fai):
        fail("FASTA index (.fai) not found: {0} "
             "(igv-reports needs it)".format(fai))

    header, rows = read_tsv(args.sites_tsv)

    title = "{0} {1}".format(args.sample, args.mode)

    if not rows:
        msg = ("No {0} sites for this sample (input had no data rows). "
               "This is an expected outcome, not an error."
               .format(args.mode))
        eprint("[info] " + msg)
        write_placeholder_html(args.out_html, title, msg)
        print("sample={0} mode={1} sites=0 output={2}"
              .format(args.sample, args.mode, args.out_html))
        return 0

    # Build the sites file in a temp location.
    tmp_dir = tempfile.mkdtemp(prefix="igv_sites_")
    sites_path = os.path.join(tmp_dir, "{0}.{1}.sites.tsv"
                              .format(args.sample, args.mode))

    if args.mode == "somatic":
        if SOMATIC_CHROM_COL not in header or SOMATIC_POS_COL not in header:
            fail("somatic TSV lacks '{0}'/'{1}' columns; got: {2}"
                 .format(SOMATIC_CHROM_COL, SOMATIC_POS_COL, header))
        _, n_sites = build_somatic_sites(rows, header, sites_path)
    else:
        _, n_sites = build_translocation_sites(rows, header, sites_path)

    if n_sites == 0:
        msg = ("Input had rows but none yielded a usable coordinate "
               "(check chrom/pos columns).")
        eprint("[warn] " + msg)
        write_placeholder_html(args.out_html, title, msg)
        print("sample={0} mode={1} sites=0 output={2}"
              .format(args.sample, args.mode, args.out_html))
        return 0

    # Defensive (re)index, then render.
    ensure_bam_index(args.bam, samtools=args.samtools)
    run_create_report(
        sites_path=sites_path,
        fasta=args.fasta,
        bam=args.bam,
        out_html=args.out_html,
        flanking=args.flanking,
        title=title,
        create_report=args.create_report,
    )

    print("sample={0} mode={1} sites={2} output={3}"
          .format(args.sample, args.mode, n_sites, args.out_html))
    return 0


if __name__ == "__main__":
    sys.exit(main())
