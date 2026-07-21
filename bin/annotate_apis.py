#!/usr/bin/env python3
"""
annotate_apis.py
---------------------------------------------------------------------------
Annotate the MM adaptive-WGS clinical tables with external knowledge bases:

  somatic (hg38)        v6_clinical.tsv  -> adds GeneBe (ACMG, gnomAD, ClinVar,
                        predictors) and OncoKB (oncogenicity, therapeutic level)
                        columns, keyed on chrom/pos/ref/alt.

  translocations (T2T)  mm_annotated.tsv -> adds OncoKB structural-variant
                        (fusion) annotation keyed on gene_a/gene_b. GeneBe does
                        not apply to SVs and is skipped.

Credentials are read from the tspipe credentials file (a Nextflow params{}
block): genebe_enabled / genebe_user / genebe_key / oncokb_enabled /
oncokb_token. Values are never printed. Environment variables GENEBE_USER,
GENEBE_KEY, ONCOKB_TOKEN override the file if set.

Failure policy (per request): on any unreachable API, non-200, or not-found,
the annotation columns are left empty and processing continues. A per-row
status column (genebe_status / oncokb_status: ok | not_found | error |
skipped | disabled) records the outcome so annotation coverage is visible
without the report breaking on a timeout.

Responses are cached on disk so reruns do not re-hit the APIs. --dry-run
prints request shapes (auth redacted) and makes no calls.

Dependencies: Python standard library only (urllib). No sample-specific
finding is hardcoded; tumor type is a study-level parameter (default OncoTree
MM). Gene/region names come from the input tables.
---------------------------------------------------------------------------
"""

import argparse
import base64
import csv
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request


GENEBE_URL = "https://api.genebe.net/cloud/api-public/v1/variants"
ONCOKB_BASE = "https://www.oncokb.org/api/v1"


def eprint(*a):
    print(*a, file=sys.stderr)


# ---------------------------------------------------------------------------
# Credentials
# ---------------------------------------------------------------------------

def parse_credentials(path):
    """Parse the tspipe params{} block (or KEY=VALUE lines). Returns a dict.
    Environment variables override file values. Values are not logged."""
    creds = {"genebe_enabled": False, "genebe_user": "", "genebe_key": "",
             "oncokb_enabled": False, "oncokb_token": ""}
    if path and os.path.isfile(path):
        with open(path) as fh:
            for raw in fh:
                line = raw.strip()
                if not line or line.startswith("//") or line.startswith("#"):
                    continue
                m = re.match(r"([A-Za-z_][A-Za-z0-9_]*)\s*[:=]\s*(.+)", line)
                if not m:
                    continue
                key, val = m.group(1), m.group(2).strip().rstrip(";").strip()
                val = val.strip("'\"")
                if key in creds:
                    if key.endswith("_enabled"):
                        creds[key] = val.lower() in ("true", "1", "yes")
                    else:
                        creds[key] = val
    # environment overrides (never logged)
    creds["genebe_user"] = os.environ.get("GENEBE_USER", creds["genebe_user"])
    creds["genebe_key"] = os.environ.get("GENEBE_KEY", creds["genebe_key"])
    creds["oncokb_token"] = os.environ.get("ONCOKB_TOKEN", creds["oncokb_token"])
    if os.environ.get("GENEBE_KEY"):
        creds["genebe_enabled"] = True
    if os.environ.get("ONCOKB_TOKEN"):
        creds["oncokb_enabled"] = True
    return creds


# ---------------------------------------------------------------------------
# HTTP (stdlib; injectable for testing)
# ---------------------------------------------------------------------------

def _http_json(method, url, headers, body=None, timeout=30, retries=3):
    """Return (status_code, parsed_json_or_None). Retries on 429/5xx with
    exponential backoff. Never raises for HTTP errors; raises only on a final
    network failure, which callers treat as 'error'."""
    data = json.dumps(body).encode() if body is not None else None
    backoff = 2
    last_exc = None
    for attempt in range(retries + 1):
        req = urllib.request.Request(url, data=data, method=method)
        for k, v in headers.items():
            req.add_header(k, v)
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                payload = resp.read().decode()
                return resp.status, (json.loads(payload) if payload else None)
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(backoff)
                backoff *= 2
                continue
            return e.code, None
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last_exc = e
            if attempt < retries:
                time.sleep(backoff)
                backoff *= 2
                continue
            raise last_exc
    return 0, None


# allow tests / callers to swap the transport
HTTP = _http_json


# ---------------------------------------------------------------------------
# Field extraction helpers (defensive: tolerate schema variation)
# ---------------------------------------------------------------------------

def first_present(d, keys):
    for k in keys:
        if isinstance(d, dict) and d.get(k) not in (None, ""):
            return d.get(k)
    return ""


GENEBE_FIELDS = {
    "genebe_acmg": ["acmg_classification"],
    "genebe_acmg_criteria": ["acmg_criteria"],
    "genebe_clinvar": ["clinvar_classification", "clinvar_disease"],
    "genebe_gnomad_af": ["gnomad_exomes_af", "gnomad_genomes_af",
                         "frequency_gnomad", "gnomad_af", "frequency"],
    "genebe_effect": ["effect", "consequences", "consequence"],
}

GENEBE_COLS = list(GENEBE_FIELDS.keys()) + ["genebe_status"]
ONCOKB_MUT_COLS = ["oncokb_oncogenic", "oncokb_effect",
                   "oncokb_highest_level", "oncokb_status"]
ONCOKB_SV_COLS = ["oncokb_sv_oncogenic", "oncokb_sv_effect",
                  "oncokb_sv_highest_level", "oncokb_status"]


def norm_chrom(c):
    return str(c)[3:] if str(c).lower().startswith("chr") else str(c)


# ---------------------------------------------------------------------------
# GeneBe (somatic)
# ---------------------------------------------------------------------------

def genebe_annotate(variants, creds, genome, cache, cfg):
    """variants: list of dict with chrom/pos/ref/alt. Returns list of dicts
    (one per variant) with GENEBE_COLS. Soft-fail to empty + status."""
    out = [dict.fromkeys(GENEBE_COLS, "") for _ in variants]
    if not creds["genebe_enabled"] or not creds["genebe_key"]:
        for r in out:
            r["genebe_status"] = "disabled"
        return out

    # split into cached vs to-query
    to_query, q_idx = [], []
    for i, v in enumerate(variants):
        key = f"genebe|{genome}|{norm_chrom(v['chrom'])}:{v['pos']}:{v['ref']}:{v['alt']}"
        if key in cache:
            _fill_genebe_row(out[i], cache[key])
        else:
            to_query.append(v)
            q_idx.append(i)

    if not to_query:
        return out

    if cfg["dry_run"]:
        eprint(f"[dry-run] GeneBe POST {GENEBE_URL}?genome={genome} "
               f"({len(to_query)} variants, Basic auth <redacted>)")
        for i in q_idx:
            out[i]["genebe_status"] = "dry_run"
        return out

    token = base64.b64encode(
        f"{creds['genebe_user']}:{creds['genebe_key']}".encode()).decode()
    headers = {"Authorization": f"Basic {token}",
               "Content-Type": "application/json", "Accept": "application/json"}

    for start in range(0, len(to_query), cfg["batch_size"]):
        chunk = to_query[start:start + cfg["batch_size"]]
        body = [{"chr": norm_chrom(v["chrom"]), "pos": int(v["pos"]),
                 "ref": v["ref"], "alt": v["alt"]} for v in chunk]
        url = f"{GENEBE_URL}?genome={genome}"
        try:
            status, obj = HTTP("POST", url, headers, body,
                               cfg["timeout"], cfg["retries"])
        except Exception:
            status, obj = 0, None
        results = (obj.get("variants") if isinstance(obj, dict) else obj) or []
        for j, v in enumerate(chunk):
            gi = q_idx[start + j]
            if status == 200 and j < len(results) and isinstance(results[j], dict):
                _fill_genebe_row(out[gi], results[j])
                key = (f"genebe|{genome}|{norm_chrom(v['chrom'])}:"
                       f"{v['pos']}:{v['ref']}:{v['alt']}")
                cache[key] = results[j]
            elif status == 200:
                out[gi]["genebe_status"] = "not_found"
            else:
                out[gi]["genebe_status"] = "error"
    return out


def _fill_genebe_row(row, resp):
    for col, keys in GENEBE_FIELDS.items():
        val = first_present(resp, keys)
        if isinstance(val, (list, dict)):
            val = json.dumps(val)
        row[col] = val
    row["genebe_status"] = "ok"


# ---------------------------------------------------------------------------
# OncoKB (somatic mutations + structural variants)
# ---------------------------------------------------------------------------

def _oncokb_headers(creds):
    return {"Authorization": f"Bearer {creds['oncokb_token']}",
            "Content-Type": "application/json", "Accept": "application/json"}


def _genomic_location(v):
    """chrom,start,end,ref,alt with chrom de-prefixed. SNV exact;
    indels best-effort (status will flag misses)."""
    chrom = norm_chrom(v["chrom"])
    pos = int(v["pos"])
    ref, alt = v["ref"], v["alt"]
    end = pos + max(0, len(ref) - 1)
    return f"{chrom},{pos},{end},{ref},{alt}"


def oncokb_annotate_mutations(variants, creds, ref_genome, tumor, cache, cfg):
    out = [dict.fromkeys(ONCOKB_MUT_COLS, "") for _ in variants]
    if not creds["oncokb_enabled"] or not creds["oncokb_token"]:
        for r in out:
            r["oncokb_status"] = "disabled"
        return out

    to_query, q_idx, locs = [], [], []
    for i, v in enumerate(variants):
        loc = _genomic_location(v)
        key = f"oncokb_mut|{ref_genome}|{tumor}|{loc}"
        if key in cache:
            _fill_oncokb_row(out[i], cache[key], ONCOKB_MUT_COLS)
        else:
            to_query.append(v); q_idx.append(i); locs.append(loc)

    if not to_query:
        return out

    if cfg["dry_run"]:
        eprint(f"[dry-run] OncoKB POST {ONCOKB_BASE}/annotate/mutations/"
               f"byGenomicChange ({len(to_query)} variants, ref={ref_genome}, "
               f"tumor={tumor}, Bearer <redacted>)")
        for i in q_idx:
            out[i]["oncokb_status"] = "dry_run"
        return out

    url = f"{ONCOKB_BASE}/annotate/mutations/byGenomicChange"
    headers = _oncokb_headers(creds)
    for start in range(0, len(to_query), cfg["batch_size"]):
        chunk_idx = q_idx[start:start + cfg["batch_size"]]
        chunk_locs = locs[start:start + cfg["batch_size"]]
        body = [{"genomicLocation": loc, "referenceGenome": ref_genome,
                 "tumorType": tumor} for loc in chunk_locs]
        try:
            status, obj = HTTP("POST", url, headers, body,
                               cfg["timeout"], cfg["retries"])
        except Exception:
            status, obj = 0, None
        results = obj if isinstance(obj, list) else []
        for j, gi in enumerate(chunk_idx):
            if status == 200 and j < len(results) and isinstance(results[j], dict):
                _fill_oncokb_row(out[gi], results[j], ONCOKB_MUT_COLS)
                cache[f"oncokb_mut|{ref_genome}|{tumor}|{chunk_locs[j]}"] = results[j]
            elif status == 200:
                out[gi]["oncokb_status"] = "not_found"
            else:
                out[gi]["oncokb_status"] = "error"
    return out


def oncokb_annotate_svs(pairs, creds, ref_genome, tumor, cache, cfg):
    """pairs: list of dict with gene_a/gene_b. Returns list with ONCOKB_SV_COLS."""
    out = [dict.fromkeys(ONCOKB_SV_COLS, "") for _ in pairs]
    if not creds["oncokb_enabled"] or not creds["oncokb_token"]:
        for r in out:
            r["oncokb_status"] = "disabled"
        return out

    def clean(g):
        g = (g or "").strip()
        return re.sub(r"_locus$", "", g, flags=re.I)

    to_query, q_idx, bodies = [], [], []
    for i, p in enumerate(pairs):
        ga, gb = clean(p.get("gene_a")), clean(p.get("gene_b"))
        if not ga or not gb or ga in (".", "NA") or gb in (".", "NA"):
            out[i]["oncokb_status"] = "skipped"
            continue
        key = f"oncokb_sv|{ref_genome}|{tumor}|{ga}__{gb}"
        if key in cache:
            _fill_oncokb_row(out[i], cache[key], ONCOKB_SV_COLS)
        else:
            to_query.append((ga, gb)); q_idx.append(i)
            bodies.append({"geneA": {"hugoSymbol": ga},
                           "geneB": {"hugoSymbol": gb},
                           "structuralVariantType": "TRANSLOCATION",
                           "functionalFusion": True,
                           "tumorType": tumor, "referenceGenome": ref_genome})

    if not bodies:
        return out

    if cfg["dry_run"]:
        eprint(f"[dry-run] OncoKB POST {ONCOKB_BASE}/annotate/structuralVariants "
               f"({len(bodies)} pairs, ref={ref_genome}, tumor={tumor}, "
               f"Bearer <redacted>)")
        for i in q_idx:
            out[i]["oncokb_status"] = "dry_run"
        return out

    url = f"{ONCOKB_BASE}/annotate/structuralVariants"
    headers = _oncokb_headers(creds)
    try:
        status, obj = HTTP("POST", url, headers, bodies,
                           cfg["timeout"], cfg["retries"])
    except Exception:
        status, obj = 0, None
    results = obj if isinstance(obj, list) else []
    for j, gi in enumerate(q_idx):
        ga, gb = to_query[j]
        if status == 200 and j < len(results) and isinstance(results[j], dict):
            _fill_oncokb_row(out[gi], results[j], ONCOKB_SV_COLS)
            cache[f"oncokb_sv|{ref_genome}|{tumor}|{ga}__{gb}"] = results[j]
        elif status == 200:
            out[gi]["oncokb_status"] = "not_found"
        else:
            out[gi]["oncokb_status"] = "error"
    return out


def _fill_oncokb_row(row, resp, cols):
    onc = resp.get("oncogenic", "")
    eff = ""
    me = resp.get("mutationEffect")
    if isinstance(me, dict):
        eff = me.get("knownEffect", "")
    level = resp.get("highestSensitiveLevel") or resp.get("highestResistanceLevel") or ""
    # cols layout: [oncogenic, effect, highest_level, status]
    row[cols[0]] = onc
    row[cols[1]] = eff
    row[cols[2]] = level
    row[cols[3]] = "ok" if (onc not in ("", "Unknown") or resp.get("variantExist")
                            or resp.get("geneExist")) else "not_found"


# ---------------------------------------------------------------------------
# TSV plumbing
# ---------------------------------------------------------------------------

def read_tsv(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    if not rows:
        return [], []
    header = rows[0]
    body = [dict(zip(header, r)) for r in rows[1:] if any(c.strip() for c in r)]
    return header, body


def write_tsv(path, header, rows, new_cols):
    out_header = header + [c for c in new_cols if c not in header]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(out_header)
        for r in rows:
            w.writerow([r.get(c, "") for c in out_header])


def load_cache(path):
    if path and os.path.isfile(path):
        try:
            with open(path) as fh:
                return json.load(fh)
        except (ValueError, OSError):
            return {}
    return {}


def save_cache(path, cache):
    if not path:
        return
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(cache, fh)


def status_tally(rows, col):
    counts = {}
    for r in rows:
        counts[r.get(col, "")] = counts.get(r.get(col, ""), 0) + 1
    return ", ".join(f"{k or '-'}={v}" for k, v in sorted(counts.items()))


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Annotate MM tables via GeneBe + OncoKB.")
    ap.add_argument("--somatic-tsv", default="", help="v6_clinical.tsv (hg38)")
    ap.add_argument("--sv-tsv", default="", help="mm_annotated.tsv (translocations)")
    ap.add_argument("--somatic-out", default="", help="default: <in>.annotated.tsv")
    ap.add_argument("--sv-out", default="", help="default: <in>.annotated.tsv")
    ap.add_argument("--credentials-file",
                    default="/home/hemat/.config/nf-core-tspipe/credentials.config")
    ap.add_argument("--genome", default="hg38", help="GeneBe genome (hg38/hg19)")
    ap.add_argument("--ref-genome", default="GRCh38", help="OncoKB referenceGenome")
    ap.add_argument("--tumor-type", default="MM", help="OncoTree code (default MM)")
    ap.add_argument("--cache", default="", help="response cache JSON (default: alongside output)")
    ap.add_argument("--batch-size", type=int, default=500)
    ap.add_argument("--timeout", type=int, default=30)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    if not args.somatic_tsv and not args.sv_tsv:
        eprint("ERROR: provide --somatic-tsv and/or --sv-tsv")
        sys.exit(2)

    creds = parse_credentials(args.credentials_file)
    eprint(f"[info] GeneBe enabled: {creds['genebe_enabled']} | "
           f"OncoKB enabled: {creds['oncokb_enabled']} "
           f"(values not shown)")

    cfg = {"dry_run": args.dry_run, "batch_size": args.batch_size,
           "timeout": args.timeout, "retries": args.retries}

    cache_path = args.cache
    if not cache_path:
        anchor = args.somatic_tsv or args.sv_tsv
        cache_path = os.path.join(os.path.dirname(os.path.abspath(anchor)),
                                  ".api_cache.json")
    cache = load_cache(cache_path)

    # --- somatic ---
    if args.somatic_tsv:
        header, rows = read_tsv(args.somatic_tsv)
        variants = [{"chrom": r.get("chrom"), "pos": r.get("pos"),
                     "ref": r.get("ref"), "alt": r.get("alt")} for r in rows]
        gb = genebe_annotate(variants, creds, args.genome, cache, cfg)
        ok = oncokb_annotate_mutations(variants, creds, args.ref_genome,
                                       args.tumor_type, cache, cfg)
        for i, r in enumerate(rows):
            r.update(gb[i])
            for c in ONCOKB_MUT_COLS:
                r[c] = ok[i][c]
        out = args.somatic_out or re.sub(r"\.tsv$", ".annotated.tsv", args.somatic_tsv)
        write_tsv(out, header, rows, GENEBE_COLS + ONCOKB_MUT_COLS)
        eprint(f"[done] somatic -> {out}")
        eprint(f"       GeneBe: {status_tally(rows, 'genebe_status')}")
        eprint(f"       OncoKB: {status_tally(rows, 'oncokb_status')}")

    # --- translocations ---
    if args.sv_tsv:
        header, rows = read_tsv(args.sv_tsv)
        pairs = [{"gene_a": r.get("gene_a"), "gene_b": r.get("gene_b")} for r in rows]
        sv = oncokb_annotate_svs(pairs, creds, args.ref_genome,
                                 args.tumor_type, cache, cfg)
        for i, r in enumerate(rows):
            for c in ONCOKB_SV_COLS:
                r[c] = sv[i][c]
        out = args.sv_out or re.sub(r"\.tsv$", ".annotated.tsv", args.sv_tsv)
        write_tsv(out, header, rows, ONCOKB_SV_COLS)
        eprint(f"[done] translocations -> {out}")
        eprint(f"       OncoKB SV: {status_tally(rows, 'oncokb_status')}")

    if not args.dry_run:
        save_cache(cache_path, cache)
        eprint(f"[info] cache: {cache_path}")


if __name__ == "__main__":
    main()
