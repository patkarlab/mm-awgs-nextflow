#!/usr/bin/env python3
"""
build_sample_report.py
---------------------------------------------------------------------------
Assemble ONE per-patient report at report/<sample>/<sample>.report.html that
ties together everything already produced by the pipeline:

  - QC panel (the sample_qc.py *.qc.html fragment) injected as create_report
    --header, per track
  - the annotated, filterable variant tables (igv-reports --tabulator) showing
    GeneBe + OncoKB columns via --info-columns
  - an Include/Exclude checkbox column + "Export included (TSV)" button added
    via create_report --footer JS (no template fork)
  - the native, working IGV pileup (the create_report page itself; never a link)

Two reference tracks become two tabs in one self-contained file:
  somatic        hg38 BAM + v6_clinical.annotated.tsv
  translocations T2T  BAM + mm_annotated.annotated.tsv (each SV row exploded
                 into its two breakends, sharing one event id)

The tabbed file embeds each track's full create_report page via lazy-loaded
iframe srcdoc (held in a text/plain element, assigned on first tab activation
so IGV initializes in a visible, correctly-sized frame). The result is one
portable HTML with no external dependencies.

stdlib only; needs create_report (igv-reports) on PATH (awgs_sv). No variant,
gene-pair, or finding is hardcoded; only column names and the panel are read.
---------------------------------------------------------------------------
"""

import argparse
import csv
import glob
import html
import os
import re
import subprocess
import sys
import tempfile


def eprint(*a):
    print(*a, file=sys.stderr)


# Columns surfaced in each track's table (only those present are used).
SOMATIC_INFO = ["gene", "csq_primary", "impact", "tumor_af_pct", "DP", "Filter",
                "genebe_acmg", "genebe_clinvar", "genebe_gnomad_af",
                "oncokb_oncogenic", "oncokb_effect", "oncokb_highest_level",
                "oncokb_status", "include"]
SV_INFO = ["event", "breakend", "sv_type", "gene_a", "gene_b",
           "known_mm_pair", "known_freq", "n_callers", "support_reads",
           "oncokb_sv_oncogenic", "oncokb_sv_effect", "oncokb_sv_highest_level",
           "oncokb_status", "include"]

FALSY = {"", ".", "0", "false", "no", "n", "na", "none", "nan"}


def truthy(v):
    return str(v).strip().lower() not in FALSY


def read_tsv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return [], []
    return rows[0], [dict(zip(rows[0], r)) for r in rows[1:] if any(c.strip() for c in r)]


def write_sites(path, columns, records):
    """columns[0:3] must be chrom/begin/end; rest are info columns."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(columns)
        for rec in records:
            w.writerow([rec.get(c, "") for c in columns])


def as_int(x):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return 0


# ---------------------------------------------------------------------------
# Sites construction
# ---------------------------------------------------------------------------

def somatic_sites(header, rows):
    """chrom/pos/pos + info columns. Pre-seed include for likely-reportable."""
    info = [c for c in SOMATIC_INFO if c in header or c == "include"]
    cols = ["chrom", "begin", "end"] + [c for c in info if c != "include"] + ["include"]
    out = []
    for r in rows:
        rec = dict(r)
        rec["begin"] = r.get("pos", "")
        rec["end"] = r.get("pos", "")
        path = (truthy(r.get("Filter", "")) and r.get("Filter", "").upper() == "PASS")
        onco = str(r.get("oncokb_oncogenic", "")).lower()
        acmg = str(r.get("genebe_acmg", "")).lower()
        reportable = path and ("oncogenic" in onco or "pathogenic" in acmg)
        rec["include"] = "1" if reportable else "0"
        out.append(rec)
    return cols, out


def sv_sites(header, rows, sv_max):
    """Explode each translocation into two breakend sites sharing an event id.
    Only interchromosomal BND/TRA events (or known MM pairs) are kept; INS/DEL/
    DUP/INV are excluded -- focal copy number is handled by depth-of-coverage on
    the panel regions, broad CNV by ichorCNA. Restrict to reportable rows and
    cap at sv_max."""
    transloc_types = {"BND", "TRA", "TRANSLOCATION", "CTX"}

    def is_translocation(r):
        t = str(r.get("sv_type", "")).upper()
        ca, cb = r.get("chrom_a", ""), r.get("chrom_b", "")
        interchrom = truthy(ca) and truthy(cb) and ca != cb
        return (t in transloc_types and interchrom) or truthy(r.get("known_mm_pair", ""))

    rows = [r for r in rows if is_translocation(r)]

    def priority(r):
        known = 0 if truthy(r.get("known_mm_pair", "")) else 1
        ok = 0 if str(r.get("oncokb_status", "")) == "ok" else 1
        return (known, ok, -as_int(r.get("n_callers")), -as_int(r.get("support_reads")))

    reportable = [r for r in rows
                  if truthy(r.get("known_mm_pair", ""))
                  or str(r.get("oncokb_status", "")) == "ok"
                  or as_int(r.get("n_callers")) >= 2]
    reportable.sort(key=priority)
    reportable = reportable[:sv_max]

    info = [c for c in SV_INFO if c in header or c in ("event", "breakend", "include")]
    cols = ["chrom", "begin", "end"] + [c for c in info if c not in ("event", "breakend", "include")] \
        + ["event", "breakend", "include"]
    out = []
    for i, r in enumerate(reportable):
        event = r.get("sv_id") or f"event{i+1}"
        inc = "1" if (truthy(r.get("known_mm_pair", ""))
                      or str(r.get("oncokb_status", "")) == "ok") else "0"
        ends = [("A", r.get("chrom_a"), r.get("pos_a"))]
        if r.get("chrom_b") and truthy(r.get("chrom_b")) and r.get("pos_b") and truthy(r.get("pos_b")):
            ends.append(("B", r.get("chrom_b"), r.get("pos_b")))
        for tag, chrom, pos in ends:
            if not chrom or not pos:
                continue
            rec = dict(r)
            rec["chrom"] = chrom
            rec["begin"] = pos
            rec["end"] = pos
            rec["event"] = event
            rec["breakend"] = tag
            rec["include"] = inc
            out.append(rec)
    return cols, out, len(reportable)


# ---------------------------------------------------------------------------
# create_report invocation
# ---------------------------------------------------------------------------

FOOTER_JS = """
<div id="curation-bar" style="position:sticky;top:0;background:#fff;
 border-bottom:1px solid #ddd;padding:6px 10px;z-index:50;
 font:13px -apple-system,Segoe UI,Roboto,Arial,sans-serif;">
  <button id="exportIncluded" type="button">Export included (TSV)</button>
  <span id="inclCount" style="margin-left:10px;color:#444;"></span>
  <span style="margin-left:10px;color:#888;">tick &ldquo;Include&rdquo; to mark a variant for the report</span>
</div>
<script>
(function(){
  var tries=0;
  function updateCount(){
    try{var inc=table.getData().filter(function(r){return r.__include;}).length;
      document.getElementById('inclCount').textContent=inc+' / '+table.getData().length+' included';}catch(e){}
  }
  var iv=setInterval(function(){
    tries++;
    if(typeof table==='undefined'||!table||!table.getColumns){ if(tries>400){clearInterval(iv);} return; }
    try{
      if(table.getData().length===0 && tries<400){ return; }
      table.getRows().forEach(function(r){var d=r.getData();
        if(d.__include===undefined){r.update({__include:String(d.include)==='1'||d.include===true});}});
      if(!table.getColumns().some(function(c){return c.getField()==='__include';})){
        table.addColumn({title:"Include",field:"__include",formatter:"tickCross",
          hozAlign:"center",headerSort:false,editor:true,width:90,frozen:true,
          cellEdited:updateCount},true);
      }
      clearInterval(iv); updateCount();
    }catch(e){ if(tries>400){clearInterval(iv);} }
  },50);
  document.addEventListener('click',function(ev){
    if(ev.target && ev.target.id==='exportIncluded'){
      if(typeof table==='undefined'||!table){return;}
      var rows=table.getData().filter(function(r){return r.__include;});
      if(!rows.length){alert('No variants ticked for inclusion.');return;}
      var cols=Object.keys(rows[0]).filter(function(k){return k!=='__include'&&k!=='unique_id';});
      var tsv=cols.join('\\t')+'\\n'+rows.map(function(r){return cols.map(function(c){
        return (r[c]===undefined||r[c]===null)?'':String(r[c]);}).join('\\t');}).join('\\n');
      var blob=new Blob([tsv],{type:'text/tab-separated-values'});
      var a=document.createElement('a'); a.href=URL.createObjectURL(blob);
      a.download='curated_included.tsv'; document.body.appendChild(a); a.click(); a.remove();
    }
  });
})();
</script>
"""


def run_create_report(sites_tsv, fasta, bam, header_html, out_html, flanking,
                      info_cols, title, create_report_exe):
    footer_path = None
    header_path = header_html
    with tempfile.NamedTemporaryFile("w", suffix=".footer.html", delete=False) as fh:
        fh.write(FOOTER_JS)
        footer_path = fh.name
    cmd = [create_report_exe, sites_tsv, fasta,
           "--sequence", "1", "--begin", "2", "--end", "3",
           "--flanking", str(flanking), "--tracks", bam,
           "--tabulator", "--info-columns"] + info_cols + \
          ["--footer", footer_path, "--title", title, "--output", out_html]
    if header_path and os.path.isfile(header_path):
        cmd += ["--header", header_path]
    eprint("[create_report] " + " ".join(cmd[:6]) + " ... --output " + out_html)
    try:
        subprocess.run(cmd, check=True)
    finally:
        if footer_path and os.path.isfile(footer_path):
            os.unlink(footer_path)
    return out_html


# ---------------------------------------------------------------------------
# Tabbed stitch (lazy srcdoc; no external files, no attribute escaping)
# ---------------------------------------------------------------------------

ENDSCRIPT_SENTINEL = "%%IGV_ENDSCRIPT%%"


def neutralize(page):
    # store the raw page inside a <script type=text/plain> holder safely by
    # replacing the closing-tag prefix with a sentinel that contains no regex
    # or escaping hazards; the loader reverses it with a plain split/join.
    return re.sub(r"</script", ENDSCRIPT_SENTINEL, page, flags=re.I)


def stitch(sample, panes, out_path):
    """panes: list of (tab_id, label, page_html). First is shown by default."""
    holders, buttons, frames = [], [], []
    for i, (tid, label, page) in enumerate(panes):
        active = " active" if i == 0 else ""
        holders.append(f'<script type="text/plain" id="src-{tid}">{neutralize(page)}</script>')
        buttons.append(f'<button class="tab{active}" data-t="{tid}">{html.escape(label)}</button>')
        frames.append(f'<iframe class="pane{active}" id="pane-{tid}" '
                      f'data-loaded="0"></iframe>')
    css = """
    body{margin:0;font:14px -apple-system,Segoe UI,Roboto,Arial,sans-serif;}
    header{padding:10px 16px;border-bottom:1px solid #ddd;}
    h1{font-size:18px;margin:0;}
    .tabs{display:flex;gap:4px;padding:8px 16px 0;border-bottom:1px solid #ddd;background:#fafafa;}
    .tab{border:1px solid #ddd;border-bottom:none;background:#eef1f3;padding:7px 16px;
      border-radius:6px 6px 0 0;cursor:pointer;font-size:13px;}
    .tab.active{background:#fff;font-weight:600;}
    iframe.pane{display:none;width:100%;height:calc(100vh - 96px);border:none;}
    iframe.pane.active{display:block;}
    """
    js = """
    function loadPane(tid){
      var f=document.getElementById('pane-'+tid);
      if(f.getAttribute('data-loaded')==='1') return;
      var raw=document.getElementById('src-'+tid).textContent;
      f.srcdoc=raw.split('%%IGV_ENDSCRIPT%%').join('</script');
      f.setAttribute('data-loaded','1');
    }
    document.querySelectorAll('.tab').forEach(function(b){
      b.addEventListener('click',function(){
        var tid=b.getAttribute('data-t');
        document.querySelectorAll('.tab').forEach(function(x){x.classList.remove('active');});
        document.querySelectorAll('iframe.pane').forEach(function(x){x.classList.remove('active');});
        b.classList.add('active');
        document.getElementById('pane-'+tid).classList.add('active');
        loadPane(tid);
      });
    });
    // load the default (first) pane
    document.addEventListener('DOMContentLoaded',function(){
      var first=document.querySelector('.tab.active');
      if(first) loadPane(first.getAttribute('data-t'));
    });
    """
    page = f"""<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(sample)} report</title><style>{css}</style></head><body>
<header><h1>{html.escape(sample)} &mdash; adaptive-WGS report</h1></header>
<div class="tabs">{''.join(buttons)}</div>
{''.join(frames)}
{''.join(holders)}
<script>{js}</script>
</body></html>"""
    with open(out_path, "w") as fh:
        fh.write(page)
    return out_path


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Build the per-patient report.")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--hg38-fasta", required=True)
    ap.add_argument("--t2t-fasta", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--flanking", type=int, default=2000)
    ap.add_argument("--sv-max", type=int, default=50)
    ap.add_argument("--create-report", default="create_report")
    # explicit overrides (else discovered from the results tree)
    ap.add_argument("--somatic-tsv", default="")
    ap.add_argument("--sv-tsv", default="")
    ap.add_argument("--hg38-bam", default="")
    ap.add_argument("--t2t-bam", default="")
    ap.add_argument("--somatic-qc", default="")
    ap.add_argument("--t2t-qc", default="")
    args = ap.parse_args()

    rd, s = os.path.abspath(args.results_dir), args.sample
    rep = os.path.join(rd, "report", s)

    def disc(p, default):
        return p if p else default

    def one(pattern):
        hits = sorted(glob.glob(pattern))
        return hits[0] if hits else ""

    somatic_tsv = disc(args.somatic_tsv, one(os.path.join(rep, "somatic", "*.v6_clinical.annotated.tsv")))
    sv_tsv = disc(args.sv_tsv, os.path.join(rep, "sv", f"{s}.mm_annotated.annotated.tsv"))
    hg38_bam = disc(args.hg38_bam, os.path.join(rd, "hg38", "bams", f"{s}.hg38.bam"))
    t2t_bam = disc(args.t2t_bam, os.path.join(rd, "t2t", "bams", f"{s}.t2t.bam"))
    somatic_qc = disc(args.somatic_qc, os.path.join(rep, "qc", f"{s}.hg38.qc.html"))
    t2t_qc = disc(args.t2t_qc, os.path.join(rep, "qc", f"{s}.t2t.qc.html"))
    out = args.out or os.path.join(rep, f"{s}.report.html")

    tmpdir = tempfile.mkdtemp(prefix="report_")
    panes = []

    # --- somatic track ---
    if os.path.isfile(somatic_tsv) and os.path.isfile(hg38_bam):
        header, rows = read_tsv(somatic_tsv)
        cols, recs = somatic_sites(header, rows)
        sites = os.path.join(tmpdir, "somatic.sites.tsv")
        write_sites(sites, cols, recs)
        info = [c for c in cols[3:]]
        page_path = os.path.join(tmpdir, "somatic.page.html")
        run_create_report(sites, args.hg38_fasta, hg38_bam,
                           somatic_qc if os.path.isfile(somatic_qc) else "",
                           page_path, args.flanking, info,
                           f"{s} somatic (hg38)", args.create_report)
        with open(page_path) as fh:
            panes.append(("somatic", "Somatic (hg38)", fh.read()))
        eprint(f"[info] somatic: {len(rows)} variants")
    else:
        eprint(f"[warn] somatic inputs missing; skipping somatic tab")

    # --- translocation track ---
    if os.path.isfile(sv_tsv) and os.path.isfile(t2t_bam):
        header, rows = read_tsv(sv_tsv)
        cols, recs, n_report = sv_sites(header, rows, args.sv_max)
        if recs:
            sites = os.path.join(tmpdir, "transloc.sites.tsv")
            write_sites(sites, cols, recs)
            info = [c for c in cols[3:]]
            page_path = os.path.join(tmpdir, "transloc.page.html")
            run_create_report(sites, args.t2t_fasta, t2t_bam,
                               t2t_qc if os.path.isfile(t2t_qc) else "",
                               page_path, args.flanking, info,
                               f"{s} translocations (T2T)", args.create_report)
            with open(page_path) as fh:
                panes.append(("transloc", "Translocations (T2T)", fh.read()))
            eprint(f"[info] translocations: {n_report} reportable of {len(rows)} "
                   f"(exploded to {len(recs)} breakend sites)")
        else:
            eprint("[warn] no reportable translocations; skipping translocation tab")
    else:
        eprint(f"[warn] translocation inputs missing; skipping translocation tab")

    if not panes:
        eprint("ERROR: no tracks could be built")
        sys.exit(2)

    os.makedirs(os.path.dirname(out), exist_ok=True)
    stitch(s, panes, out)
    eprint(f"[done] report: {out}")
    print(out)


if __name__ == "__main__":
    main()
