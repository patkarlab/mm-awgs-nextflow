#!/usr/bin/env python3
"""
add_hgvs_v2.py -- add HGVSc/HGVSp columns to Clair3 VEP output.
Anchors on short escaping-free fragments to avoid backslash-layer issues.
Idempotent (.bak, refuses double-apply). See prior notes for rationale:
norm -m- before VEP fixes the multi-allelic --hgvs crash; columns are
end-appended so nothing renumbers.
"""
import argparse, datetime, os, shutil, sys

SENTINEL = "hgvs-annotation-applied-v2"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", required=True)
    args = ap.parse_args()
    repo = args.repo

    nf = os.path.join(repo, "modules/local/vep_annotate_clair3.nf")
    filt = os.path.join(repo, "bin/filter_v6_somatic_candidates.py")
    for p in (nf, filt):
        if not os.path.isfile(p):
            sys.exit(f"ERROR: not found: {p}")

    t = open(nf).read()
    if SENTINEL in t:
        print(f"SKIP (already applied): {nf}")
    else:
        orig = t
        problems = []

        # --- Edit 1: norm step. Anchor: the FIRST `tabix -p vcf -f "$pass_vcf"`.
        # We locate the exact substring regardless of quoting by searching for
        # the pass_vcf tabix; the file uses "$pass_vcf" (Groovy passes \$ -> $;
        # on disk the bytes are backslash-dollar). Match on 'pass_vcf"' tail.
        a1 = 'tabix -p vcf -f "$pass_vcf"'
        a1b = 'tabix -p vcf -f "\\$pass_vcf"'   # escaped-dollar on-disk form
        hit1 = a1 if a1 in t else (a1b if a1b in t else None)
        if hit1:
            dollar = "\\$" if hit1 == a1b else "$"
            ins = (hit1 + "\n"
                   "    norm_vcf=vep_out/${meta.id}.norm.vcf.gz\n"
                   f'    bcftools norm -m- -f ${{params.hg38_fasta}} "{dollar}pass_vcf" -O z -o "{dollar}norm_vcf"\n'
                   f'    tabix -p vcf -f "{dollar}norm_vcf"')
            t = t.replace(hit1, ins, 1)
        else:
            problems.append("norm anchor (pass_vcf tabix)")

        # --- Edit 2: VEP input -> norm. Match 'input_file "$pass_vcf"' either form.
        for old, new in [('--input_file "$pass_vcf"', '--input_file "$norm_vcf"'),
                         ('--input_file "\\$pass_vcf"', '--input_file "\\$norm_vcf"')]:
            if old in t:
                t = t.replace(old, new, 1); break
        else:
            problems.append("VEP --input_file")

        # --- Edit 3: --hgvs after --use_given_ref (escaping-free anchor).
        if "--use_given_ref" in t:
            # insert '--hgvs \' line after the use_given_ref line
            lines = t.split("\n")
            for i, ln in enumerate(lines):
                if "--use_given_ref" in ln:
                    indent = ln[:len(ln)-len(ln.lstrip())]
                    cont = " \\" if ln.rstrip().endswith("\\") else ""
                    lines.insert(i+1, f"{indent}--hgvs{cont}")
                    break
            t = "\n".join(lines)
        else:
            problems.append("--use_given_ref")

        # --- Edit 4: split-vep -f: insert %HGVSc %HGVSp before the DP token's close.
        # Escaping-free anchor: the literal '[%DP]' substring. Insert HGVS tokens
        # after it, mirroring whatever tab-escape precedes it by capturing it.
        import re
        # find '[%DP]' then the run of backslashes+t or t before %HGVS should
        # match the sep used before [%AF]. Capture the separator preceding [%AD].
        m = re.search(r'(\\+t|\bt)?\[%AD\]', t)
        sep = None
        if m and m.group(1):
            sep = m.group(1)          # e.g. '\\t' as on disk
        if "[%DP]" in t and sep:
            t = t.replace("[%DP]", f"[%DP]{sep}%HGVSc{sep}%HGVSp", 1)
        else:
            problems.append("split-vep [%DP]/sep")

        # --- Edit 5: awk print append. Escaping-free anchor: 'rc, ac, dp' at line end.
        # Append ', $21, $22' (or escaped) matching the file's $ form.
        for old, new in [("rc, ac, dp\n", "rc, ac, dp, $21, $22\n"),
                         ("rc, ac, dp\n", "rc, ac, dp, \\$21, \\$22\n")]:
            # prefer whichever $ style the file uses elsewhere
            pass
        if "rc, ac, dp" in t:
            dollar = "\\$" if '"\\$norm_vcf"' in t else "$"
            t = t.replace("rc, ac, dp\n", f"rc, ac, dp, {dollar}21, {dollar}22\n", 1)
        else:
            problems.append("awk 'rc, ac, dp'")

        # --- Edit 6: both headers. Escaping-free anchor: 'ALT_COUNT' ... 'DP' then sep+n.
        # The header ends '...ALT_COUNT<sep>DP<sep>n'. Insert hgvsc/hgvsp before the
        # final <sep>n. Capture the sep before ALT_COUNT.
        m2 = re.search(r'(\\+t|\bt)ALT_COUNT', t)
        hsep = m2.group(1) if m2 else sep
        m3 = re.search(r'(\\+n)', t)  # newline escape form
        nl = None
        # find 'DP' immediately followed by sep-less newline token in header
        if hsep:
            # header tail as bytes: DP<hsep>n  (n preceded by backslashes)
            # Build the DP...\n pattern using hsep's backslash count for \n too.
            bs = hsep.replace("t","")  # the backslashes part
            dp_nl = f"DP{bs}n"
            dp_nl_new = f"DP{hsep}hgvsc{hsep}hgvsp{bs}n"
            cnt = t.count(dp_nl)
            if cnt >= 1:
                t = t.replace(dp_nl, dp_nl_new)   # ALL (both headers)
            else:
                problems.append(f"header 'DP..n' (found {cnt})")
        else:
            problems.append("header sep")

        if problems:
            print("ABORT nf (no changes written). Unmatched anchors:")
            for p in problems: print("   -", p)
        else:
            bak = nf + ".bak_hgvs2_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            shutil.copy2(nf, bak)
            # neutralize old workaround comment + add sentinel
            t = t.replace("NO --hgvs / --hgvsg (multi-allelic crash workaround per step 12b)",
                          f"--hgvs enabled; multi-allelics split by bcftools norm -m- above  [{SENTINEL}]", 1)
            open(nf,"w").write(t)
            print(f"PATCHED: {nf}\n  backup: {bak}")

    # --- Edit 7: filter PREFERRED_COLS ---
    tf = open(filt).read()
    if '"hgvsc"' in tf and '"hgvsp"' in tf:
        print(f"SKIP (hgvs cols already present): {filt}")
    else:
        a = '"gene", "panel_label", "transcript", "biotype", "canonical",'
        if a in tf:
            bak = filt + ".bak_hgvs2_" + datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            shutil.copy2(filt, bak)
            tf = tf.replace(a, '"gene", "panel_label", "transcript", "hgvsc", "hgvsp", "biotype", "canonical",', 1)
            open(filt,"w").write(tf)
            print(f"PATCHED: {filt}\n  backup: {bak}")
        else:
            print(f"NOTE: PREFERRED_COLS anchor not found in {filt}; may already be patched.")

if __name__ == "__main__":
    main()
