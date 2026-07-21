#!/usr/bin/env python3
"""
build_cohort_dashboard.py
---------------------------------------------------------------------------
Assemble a portable cohort overview dashboard for the MM adaptive-WGS
pipeline. It does NOT recompute anything: it discovers per-sample artifacts
already produced by the pipeline and ties them together into one self-
contained folder that can be moved off the server as a unit.

For each sample it reads (from the results tree):
  - ichorCNA tumor fraction + ploidy   report/<sample>/ichorcna/<sample>.params.txt
  - ichorCNA solution plot             report/<sample>/ichorcna/<sample>_genomeWide_all_sols.pdf
  - MM-annotated SVs / translocations  report/<sample>/sv/<sample>.mm_annotated.tsv
  - clinical somatic SNV/indels        report/<sample>/somatic/*.v6_clinical.tsv
  - somatic IGV view (hg38)            hg38/igv/<sample>/<sample>.somatic.html
  - translocation IGV view (T2T)       t2t/igv/<sample>/<sample>.translocations.html

Output (bundle mode, default): a 'dashboard/' folder containing
cohort_dashboard.html plus copies of each sample's IGV HTMLs and featured
ichorCNA PDF, linked by relative path so the folder is portable.

Design constraints (project):
  - Python standard library only (runs in awgs_sv; no pandas).
  - No variant, gene-pair, FISH finding, or expected karyotype is hardcoded.
    The script renders whatever rows the upstream filters produced; only
    column NAMES are referenced.
  - Samples are referenced by sequencing ID only (discovered from the tree).
---------------------------------------------------------------------------
"""

import argparse
import csv
import glob
import html
import os
import shutil
import sys


# --- Columns surfaced from the clinical somatic TSV (by name, not position) -
SOMATIC_COLS = [
    ("gene", "Gene"),
    ("csq_primary", "Consequence"),
    ("impact", "Impact"),
    ("tumor_af_pct", "Tumor AF %"),
    ("DP", "Depth"),
    ("clinvar_sig", "ClinVar"),
    ("pop_af_max", "Pop AF"),
    ("Filter", "Filter"),
]

# --- Columns surfaced from the MM-annotated SV TSV --------------------------
SV_COLS = [
    ("sv_type", "Type"),
    ("gene_a", "Gene A"),
    ("chrom_a", "Chr A"),
    ("pos_a", "Pos A"),
    ("gene_b", "Gene B"),
    ("chrom_b", "Chr B"),
    ("pos_b", "Pos B"),
    ("known_mm_pair", "Known MM pair"),
    ("known_freq", "Freq"),
    ("n_callers", "Callers"),
    ("support_reads", "Support"),
]

FALSY = {"", ".", "0", "false", "no", "n", "na", "none", "nan"}


def eprint(*a):
    print(*a, file=sys.stderr)


def is_truthy(value):
    return str(value).strip().lower() not in FALSY


def read_tsv(path):
    """Return (header_list, list_of_row_dicts). Empty file -> ([], [])."""
    if not path or not os.path.isfile(path):
        return [], []
    with open(path, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        rows = list(reader)
    if not rows:
        return [], []
    header = rows[0]
    body = [dict(zip(header, r)) for r in rows[1:] if any(c.strip() for c in r)]
    return header, body


def parse_ichor_params(path):
    """Return (tumor_fraction_float_or_None, ploidy_str). Robust to the
    header / data / trailing-sample-id layout ichorCNA writes."""
    if not os.path.isfile(path):
        return None, None
    tf, ploidy = None, None
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                # data row: <sample> <tumor fraction> <ploidy>
                try:
                    tf = float(parts[1])
                    ploidy = parts[2].strip()
                    break
                except ValueError:
                    continue
    return tf, ploidy


def first_glob(pattern):
    hits = sorted(glob.glob(pattern))
    return hits[0] if hits else None


def discover_samples(report_dir):
    """Sample IDs are the sub-directory names under report/ that carry a
    somatic/ folder (our anchor artifact)."""
    samples = []
    for entry in sorted(os.listdir(report_dir)):
        sdir = os.path.join(report_dir, entry)
        if os.path.isdir(sdir) and os.path.isdir(os.path.join(sdir, "somatic")):
            samples.append(entry)
    return samples


def rel(target, start_dir):
    """Relative path from start_dir to target (for portable links)."""
    return os.path.relpath(target, start_dir)


def copy_into(src, dest_dir, dest_name):
    """Copy src into dest_dir/dest_name; return the destination path, or None."""
    if not src or not os.path.isfile(src):
        return None
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, dest_name)
    shutil.copy2(src, dest)
    return dest


# ---------------------------------------------------------------------------
# HTML rendering (no external assets; minimal inline CSS)
# ---------------------------------------------------------------------------

CSS = """
:root { --bg:#fff; --fg:#1a1a1a; --muted:#666; --line:#e3e3e3;
        --accent:#244; --warn:#9a4a00; --chip:#f2f4f5; }
* { box-sizing: border-box; }
body { font: 14px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
       color: var(--fg); background: var(--bg); margin: 0; padding: 24px; }
h1 { font-size: 20px; margin: 0 0 4px; }
h2 { font-size: 16px; margin: 28px 0 8px; border-bottom: 2px solid var(--line);
     padding-bottom: 4px; }
h3 { font-size: 14px; margin: 16px 0 6px; color: var(--accent); }
.sub { color: var(--muted); margin: 0 0 16px; }
table { border-collapse: collapse; width: 100%; margin: 6px 0 12px; font-size: 13px; }
th, td { text-align: left; padding: 5px 8px; border-bottom: 1px solid var(--line);
         vertical-align: top; }
th { background: var(--chip); font-weight: 600; }
td.num { text-align: right; font-variant-numeric: tabular-nums; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 10px;
         background: var(--chip); font-size: 12px; margin-right: 6px; }
.badge.warn { background: #fdeede; color: var(--warn); }
.known { font-weight: 600; }
.card { border: 1px solid var(--line); border-radius: 8px; padding: 16px;
        margin: 16px 0; }
.links a { margin-right: 14px; }
details { margin: 8px 0; }
summary { cursor: pointer; font-weight: 600; color: var(--accent); }
iframe { width: 100%; height: 620px; border: 1px solid var(--line);
         border-radius: 6px; margin-top: 8px; background:#fafafa; }
.note { color: var(--muted); font-size: 12px; margin: 4px 0 0; }
.empty { color: var(--muted); font-style: italic; }
footer { color: var(--muted); font-size: 12px; margin-top: 32px;
         border-top: 1px solid var(--line); padding-top: 8px; }
"""


def esc(x):
    return html.escape("" if x is None else str(x))


def tf_display(tf):
    if tf is None:
        return '<span class="badge warn">not estimated &mdash; review</span>'
    if tf == 0:
        return '<span class="badge warn">TF 0 &mdash; non-diploid solution, review</span>'
    return f'<span class="badge">TF {tf*100:.1f}%</span>'


def render_table(header, rows, cols, limit=None, sort_key=None):
    present = [(k, label) for k, label in cols if k in header]
    if not present:
        return '<p class="empty">no columns matched.</p>'
    data = rows
    if sort_key:
        data = sorted(rows, key=sort_key)
    truncated = False
    if limit and len(data) > limit:
        data = data[:limit]
        truncated = True
    out = ["<table><thead><tr>"]
    out += [f"<th>{esc(label)}</th>" for _, label in present]
    out.append("</tr></thead><tbody>")
    for r in data:
        out.append("<tr>")
        for k, _ in present:
            v = r.get(k, "")
            cls = ' class="known"' if k == "known_mm_pair" and is_truthy(v) else ""
            out.append(f"<td{cls}>{esc(v)}</td>")
        out.append("</tr>")
    out.append("</tbody></table>")
    if truncated:
        out.append(f'<p class="note">showing first {limit} of {len(rows)} rows; '
                   f'full table in the SV/somatic TSV.</p>')
    return "".join(out)


def sv_sort_key(r):
    # known MM pairs first, then more callers, then more support
    known = 0 if is_truthy(r.get("known_mm_pair", "")) else 1

    def as_int(x):
        try:
            return int(float(x))
        except (ValueError, TypeError):
            return 0
    return (known, -as_int(r.get("n_callers", 0)), -as_int(r.get("support_reads", 0)))


def main():
    ap = argparse.ArgumentParser(description="Build the MM aWGS cohort dashboard.")
    ap.add_argument("--results-dir", required=True,
                    help="Pipeline results dir (contains report/, hg38/, t2t/)")
    ap.add_argument("--out-dir", default=None,
                    help="Dashboard output dir (default: <results>/report/dashboard)")
    ap.add_argument("--title", default="MM adaptive-WGS cohort dashboard")
    ap.add_argument("--sv-preview", type=int, default=15,
                    help="Max SV rows previewed per sample (default 15)")
    args = ap.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    report_dir = os.path.join(results_dir, "report")
    if not os.path.isdir(report_dir):
        eprint(f"ERROR: no report/ under {results_dir}")
        sys.exit(2)

    out_dir = os.path.abspath(args.out_dir) if args.out_dir \
        else os.path.join(report_dir, "dashboard")
    igv_out = os.path.join(out_dir, "igv")
    pdf_out = os.path.join(out_dir, "ichorcna")
    os.makedirs(out_dir, exist_ok=True)

    samples = discover_samples(report_dir)
    if not samples:
        eprint(f"ERROR: no samples (no */somatic/ under {report_dir})")
        sys.exit(2)
    eprint(f"[info] {len(samples)} sample(s): {', '.join(samples)}")

    cohort_rows = []     # for the overview table
    sample_blocks = []   # per-sample HTML

    for s in samples:
        sdir = os.path.join(report_dir, s)

        tf, ploidy = parse_ichor_params(
            os.path.join(sdir, "ichorcna", f"{s}.params.txt"))

        sv_path = first_glob(os.path.join(sdir, "sv", "*.mm_annotated.tsv"))
        sv_header, sv_rows = read_tsv(sv_path)
        n_sv = len(sv_rows)
        n_known = sum(1 for r in sv_rows if is_truthy(r.get("known_mm_pair", "")))

        som_path = first_glob(os.path.join(sdir, "somatic", "*.v6_clinical.tsv"))
        som_header, som_rows = read_tsv(som_path)
        n_mut = len(som_rows)

        # copy portable artifacts into the bundle, link by relative path
        pdf_src = first_glob(os.path.join(sdir, "ichorcna",
                                          "*_genomeWide_all_sols.pdf"))
        pdf_dest = copy_into(pdf_src, pdf_out, f"{s}.genomeWide_all_sols.pdf")

        som_html_src = os.path.join(results_dir, "hg38", "igv", s, f"{s}.somatic.html")
        trl_html_src = os.path.join(results_dir, "t2t", "igv", s, f"{s}.translocations.html")
        som_html_dest = copy_into(som_html_src, igv_out, f"{s}.somatic.html")
        trl_html_dest = copy_into(trl_html_src, igv_out, f"{s}.translocations.html")

        cohort_rows.append({
            "sample": s, "tf": tf, "ploidy": ploidy,
            "n_sv": n_sv, "n_known": n_known, "n_mut": n_mut,
        })

        # ---- per-sample block ----
        b = [f'<div class="card" id="{esc(s)}">']
        b.append(f"<h3>{esc(s)}</h3>")
        b.append("<p>")
        b.append(tf_display(tf))
        b.append(f'<span class="badge">ploidy {esc(ploidy) if ploidy else "NA"}</span>')
        b.append(f'<span class="badge">{n_known} known pair(s)</span>')
        b.append(f'<span class="badge">{n_sv} SV / {n_mut} clinical mut</span>')
        b.append("</p>")

        # links to the interactive IGV views
        b.append('<p class="links">')
        if trl_html_dest:
            b.append(f'<a href="{esc(rel(trl_html_dest, out_dir))}" '
                     f'target="_blank">Translocation IGV (T2T) &#8599;</a>')
        if som_html_dest:
            b.append(f'<a href="{esc(rel(som_html_dest, out_dir))}" '
                     f'target="_blank">Somatic IGV (hg38) &#8599;</a>')
        b.append("</p>")

        # ichorCNA plot, embedded inline by relative path
        if pdf_dest:
            b.append("<details><summary>ichorCNA solution plot</summary>")
            b.append(f'<iframe src="{esc(rel(pdf_dest, out_dir))}" '
                     f'title="ichorCNA {esc(s)}"></iframe></details>')
        else:
            b.append('<p class="empty">ichorCNA plot not found.</p>')

        # SV preview
        b.append("<h4>Structural variants / translocations</h4>")
        if n_sv:
            b.append(render_table(sv_header, sv_rows, SV_COLS,
                                  limit=args.sv_preview, sort_key=sv_sort_key))
        else:
            b.append('<p class="empty">no annotated SVs.</p>')

        # somatic mutation preview (these are already PASS in-panel; show all)
        b.append("<h4>Clinical somatic SNV / indels</h4>")
        if n_mut:
            b.append(render_table(som_header, som_rows, SOMATIC_COLS))
        else:
            b.append('<p class="empty">no clinical somatic variants.</p>')

        # inline interactive IGV (expandable; same files, previewed in place)
        if trl_html_dest:
            b.append("<details><summary>Translocation IGV (inline)</summary>")
            b.append(f'<iframe src="{esc(rel(trl_html_dest, out_dir))}"></iframe>')
            b.append("</details>")
        if som_html_dest:
            b.append("<details><summary>Somatic IGV (inline)</summary>")
            b.append(f'<iframe src="{esc(rel(som_html_dest, out_dir))}"></iframe>')
            b.append("</details>")

        b.append("</div>")
        sample_blocks.append("".join(b))

    # ---- cohort overview table ----
    ov = ['<table><thead><tr>'
          '<th>Sample</th><th>Tumor fraction</th><th>Ploidy</th>'
          '<th>SVs</th><th>Known MM pairs</th><th>Clinical mutations</th>'
          '</tr></thead><tbody>']
    for r in cohort_rows:
        if r["tf"] is None:
            tf_cell = '<span class="badge warn">review</span>'
        elif r["tf"] == 0:
            tf_cell = '<span class="badge warn">0 &mdash; review</span>'
        else:
            tf_cell = f'{r["tf"]*100:.1f}%'
        ov.append(
            f'<tr><td><a href="#{esc(r["sample"])}">{esc(r["sample"])}</a></td>'
            f'<td>{tf_cell}</td><td class="num">{esc(r["ploidy"])}</td>'
            f'<td class="num">{r["n_sv"]}</td><td class="num">{r["n_known"]}</td>'
            f'<td class="num">{r["n_mut"]}</td></tr>'
        )
    ov.append("</tbody></table>")

    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{esc(args.title)}</title><style>{CSS}</style></head>
<body>
<h1>{esc(args.title)}</h1>
<p class="sub">{len(samples)} sample(s) &middot; source: {esc(results_dir)}</p>
<h2>Cohort overview</h2>
<p class="note">Tumor fraction and ploidy are the ichorCNA selected-solution
values; TF reported as 0 indicates a non-diploid solution that needs manual
review of the plot, not a 0% tumour. Genome-wide copy-number-changing events
are covered by ichorCNA; CN-LOH is out of scope for this assay.</p>
{''.join(ov)}
<h2>Per-sample detail</h2>
{''.join(sample_blocks)}
<footer>Generated by build_cohort_dashboard.py. Portable: move the whole
'{os.path.basename(out_dir)}/' folder as a unit; all links are relative.</footer>
</body></html>
"""

    out_html = os.path.join(out_dir, "cohort_dashboard.html")
    with open(out_html, "w") as fh:
        fh.write(page)
    eprint(f"[done] dashboard: {out_html}")
    print(out_html)


if __name__ == "__main__":
    main()
