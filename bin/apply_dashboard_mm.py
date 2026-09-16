#!/usr/bin/env python3
"""
apply_dashboard_mm.py  (v3, 2026-09-16 evening)

Applies the dashboard revision (dashboard_mm_v3) from dashboard_mm_v3_files.tar.gz
beside this script. v3 over v2: the bundle also collects genebaf.sites.tsv (the allele-fraction
dots were missing without it) and the arm-test line is reworded. v2 added,
over the morning's v1:

  - per-window depth + allele-fraction plot on the copy-number tab, revealed by
    clicking an assessed window's row (inline SVG from genebaf bins and sites)
  - panel window names on the per-chromosome targets track (bin/mmplot.py,
    target_labels_v1) - ICHORCNA_PLOT re-runs for this
  - a details row under every partner locus (click the row): per-caller reads
    and filters, nanomonsv, dictionary fields, flag text, breakend coordinates
  - the "Other junctions" section shows its rows in the Reportable view
  - HGVSp and other long values wrap inside the variant card

Whole-file replacement with .bak. Each replaced file must currently be the
bc3d2da, v1 or v2 version; anything else is
reported and nothing is written.

  python3 bin/apply_dashboard_mm.py --check | (apply) | --revert
"""
import argparse, hashlib, json, os, shutil, sys, tarfile
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
TARBALL = os.path.join(HERE, "dashboard_mm_v3_files.tar.gz")
M = json.loads('{"allowed": {"bin/build_report_bundle.sh": ["d995175c0b0fd3ee47d3173cc902196e2da973d9c15e422de3b64bc415b92aff", "f099918e7baf39da4667094d9a8b22216ac0e1466ceec8b5b776baf80796bcf3"], "bin/dashboard_builder/assets/css/dashboard.css": ["1b7fe49a159e3fe00797fc6552226794c639018879c49c3dfb9b16412b0192c6", "41483e9a11d27e1ab0dcd4316a7bbb28ea09243b98f6a6882c3d0dac8c8cbe96", "482a31d622288ff765038777911fac9f4e6c1256b11c6d123350f409c6e88a38"], "bin/dashboard_builder/build.py": ["31cb24804626aee5b0e7365db1d75283e77755e6b837e00a05f07aaff7801e5d", "c9c0ba0cdbcbeb05db549155647709befa3b8891d0055201e576333e55ceddc9"], "bin/dashboard_builder/templates/cohort_index.html.j2": ["1159a14b79dfcfcb9132433d5102736bae27cbe4ddc44959585104f2b7b6db0b", "fd443a8c4784ddbcc7cdf282bb6dfccf1a37984718ff4b2f3f02a6a840cd043d"], "bin/dashboard_builder/templates/copy_number_tab.html.j2": ["2eff734e735f14b867b6e3ce0b7332c7dda5ec9e28b855e773926f4c1e9f6022", "cda393cd7c27d1640c058f23ac7f5ac798fcd401e7d84d7714a7a5973a21c425", "f6fbc778914d8bbded896c7873d39f4f91226f3ee02d2658011e379ce9e47b38"], "bin/dashboard_builder/templates/sample_report.html.j2": ["4671425ce2a7aa8f0a6288191ffd000b4d132a90b82091afa51e795d736b3551", "ebce012714034f611284f16f585764a19596e30ef8dc0bf4463307a0823b8d55"], "bin/merge_translocations.py": ["7356b8d7db9def751a029805376e2b05131f2d8f18cc7f981562c52db7f29f2f", "aa4d9a36c683da3c0e35f1db438cd401b9ec386c7e90e3c38b8ab05fcbc0940c"], "bin/mmplot.py": ["8c882008165803333242954c910c417de6c2ea4a3af08e00bc62729ce64e174b", "f4f7c34ce7536c7bcc5a173b8a0411e45b4ed5f52aec2814d399585af8a6474c"], "bin/dashboard_builder/parsers/allelic.py": ["1115fe73809d37a2907158a253ed85cbef9d4f6d4f0c548060533b8405f2fb6d", "8490a256fcf2d885aee1dda5bf6ed3a80d5074cb3ae62f0cb7904e5fcc1d553e"], "bin/dashboard_builder/parsers/rearrangements.py": ["1d2e145825adbe22a14ea6256bb7ecd39ef7ac87996aa4ae30abd1f226d8f3a8", "71edc4e8e84e1aa16091bbfdd91d309cea5c78520b3aeb9915e9231b9cd04a3c"], "bin/dashboard_builder/templates/rearrangements_tab.html.j2": ["07bbbb90d08fb8471050caa2b2858a25db04eeb030303de38ddbf46c218ba594", "3cd0dac78f5bd78a56a82a116fb9120a5ab46d008e5b4a8036303b082276b512"]}, "target": {"bin/build_report_bundle.sh": "0249ff495b9b5451df01035dc8ab490870636cac2da73d270c78f0bbb28b6aa0", "bin/dashboard_builder/assets/css/dashboard.css": "1b7fe49a159e3fe00797fc6552226794c639018879c49c3dfb9b16412b0192c6", "bin/dashboard_builder/build.py": "c9c0ba0cdbcbeb05db549155647709befa3b8891d0055201e576333e55ceddc9", "bin/dashboard_builder/templates/cohort_index.html.j2": "fd443a8c4784ddbcc7cdf282bb6dfccf1a37984718ff4b2f3f02a6a840cd043d", "bin/dashboard_builder/templates/copy_number_tab.html.j2": "f66c991f630c89feaf42e3620d48589e638e289d005fb6153ba7fbf269a0ed15", "bin/dashboard_builder/templates/sample_report.html.j2": "4671425ce2a7aa8f0a6288191ffd000b4d132a90b82091afa51e795d736b3551", "bin/merge_translocations.py": "7356b8d7db9def751a029805376e2b05131f2d8f18cc7f981562c52db7f29f2f", "bin/mmplot.py": "f4f7c34ce7536c7bcc5a173b8a0411e45b4ed5f52aec2814d399585af8a6474c", "bin/dashboard_builder/parsers/allelic.py": "8490a256fcf2d885aee1dda5bf6ed3a80d5074cb3ae62f0cb7904e5fcc1d553e", "bin/dashboard_builder/parsers/rearrangements.py": "1d2e145825adbe22a14ea6256bb7ecd39ef7ac87996aa4ae30abd1f226d8f3a8", "bin/dashboard_builder/templates/rearrangements_tab.html.j2": "3cd0dac78f5bd78a56a82a116fb9120a5ab46d008e5b4a8036303b082276b512"}, "replaced": ["bin/build_report_bundle.sh", "bin/dashboard_builder/assets/css/dashboard.css", "bin/dashboard_builder/build.py", "bin/dashboard_builder/templates/cohort_index.html.j2", "bin/dashboard_builder/templates/copy_number_tab.html.j2", "bin/dashboard_builder/templates/sample_report.html.j2", "bin/merge_translocations.py", "bin/mmplot.py"], "new": ["bin/dashboard_builder/parsers/allelic.py", "bin/dashboard_builder/parsers/rearrangements.py", "bin/dashboard_builder/templates/rearrangements_tab.html.j2"]}')

def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()

def check(verbose=True):
    problems = 0
    if not os.path.isfile(TARBALL):
        print("MISSING  dashboard_mm_v3_files.tar.gz (copy it beside this script)"); return 1
    done = 0
    for rel in M["replaced"] + M["new"]:
        p = os.path.join(REPO, rel)
        if not os.path.isfile(p):
            if rel in M["new"]:
                if verbose: print(f"NEW      {rel}")
            else:
                print(f"MISSING  {rel}"); problems += 1
            continue
        h = sha(p)
        if h == M["target"][rel]:
            done += 1
            if verbose: print(f"APPLIED  {rel}")
        elif h in M["allowed"][rel]:
            if verbose: print(f"READY    {rel}")
        else:
            print(f"MODIFIED {rel}: not the bc3d2da or v1 version; tell me before applying"); problems += 1
    return 1 if problems else 0

def apply():
    if check(verbose=False) != 0:
        sys.exit("ERROR: preconditions failed; nothing changed.")
    for rel in M["replaced"]:
        p = os.path.join(REPO, rel)
        if os.path.isfile(p) and sha(p) != M["target"][rel]:
            shutil.copy2(p, p + ".bak")
    with tarfile.open(TARBALL) as tf:
        tf.extractall(REPO)
    for rel in M["replaced"] + M["new"]:
        print(f"written  {rel}")
    print("\nNext: resume the session. ICHORCNA_PLOT (labels), MERGE_TRANSLOCATIONS, IGV pages, bundle, dashboard re-run; nothing else.")

def revert():
    for rel in M["replaced"]:
        p = os.path.join(REPO, rel)
        if os.path.isfile(p + ".bak"):
            shutil.move(p + ".bak", p); print(f"restored {rel}")
    for rel in M["new"]:
        p = os.path.join(REPO, rel)
        if os.path.isfile(p):
            os.remove(p); print(f"removed  {rel}")

def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--check", action="store_true"); g.add_argument("--revert", action="store_true")
    a = ap.parse_args()
    if a.check: sys.exit(check())
    if a.revert: revert(); return
    apply()

if __name__ == "__main__":
    main()
