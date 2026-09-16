#!/usr/bin/env python3
"""
apply_sites_all.py  (sites_all_v1)

genebaf.py wrote only two-allele sites to <sample>.genebaf.sites.tsv, so on a
window with a deletion or LOH the allele-fraction plot went empty rather than
showing ~600 sites piled at 0 and 1. Every usable site is now written, with a
`call` column (het / hom_ref / hom_alt). The dashboard already colours by
allele fraction, so nothing else changes. GENEBAF re-runs.

Idempotent; --check; --revert.
"""
import argparse, os, shutil, sys
SENTINEL = "sites_all_v1"
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TARGET = "bin/genebaf.py"
EDITS = [
    ('''            if r >= args.min_alt and a >= args.min_alt:
                n_het += 1
                frac = a / n
                afs.append(frac)
                devs.append(abs(frac - 0.5))
                per_site_rows.append((name, chrom, pos, r, a, round(frac, 3)))
            else:
                n_hom += 1
''', '''            # sites_all_v1: every usable site is written, homozygous ones
            # included. A deletion or LOH empties the heterozygous band into
            # the 0/1 rows, and a reader needs to see those rows to know the
            # window was measured rather than empty.
            frac = a / n
            if r >= args.min_alt and a >= args.min_alt:
                n_het += 1
                afs.append(frac)
                devs.append(abs(frac - 0.5))
                per_site_rows.append((name, chrom, pos, r, a, round(frac, 3), "het"))
            else:
                n_hom += 1
                per_site_rows.append((name, chrom, pos, r, a, round(frac, 3),
                                      "hom_alt" if a > r else "hom_ref"))
'''),
    ('''        handle.write("gene\\tchrom\\tpos\\tref_reads\\talt_reads\\talt_fraction\\n")''',
     '''        handle.write("gene\\tchrom\\tpos\\tref_reads\\talt_reads\\talt_fraction\\tcall\\n")'''),
]
def main():
    ap = argparse.ArgumentParser(); g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true"); g.add_argument("--revert", action="store_true")
    a = ap.parse_args(); path = os.path.join(REPO, TARGET)
    if a.revert:
        if os.path.isfile(path + ".bak"): shutil.move(path + ".bak", path); print("restored", TARGET)
        return
    t = open(path, encoding="utf-8").read()
    if SENTINEL in t: print("APPLIED ", TARGET); return
    bad = [o for o, _ in EDITS if t.count(o) != 1]
    if bad: print("NO-MATCH", TARGET); sys.exit(1)
    if a.check: print("READY   ", TARGET); return
    shutil.copy2(path, path + ".bak")
    for o, n in EDITS: t = t.replace(o, n, 1)
    open(path, "w", encoding="utf-8").write(t); print("patched ", TARGET)
if __name__ == "__main__": main()
