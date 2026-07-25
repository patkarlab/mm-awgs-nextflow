# Nextflow wiring and publication

`REPORT_BUNDLE` and `DASHBOARD` are already wired in the workflow behind an
ordering gate. Three stages are added around them:

```
IGV_SNAPSHOTS -> [existing gate] -> REPORT_BUNDLE -> DASHBOARD -> EMBED_REPORT_ASSETS -> REPORT_ZIP
```

No existing call is replaced. An earlier draft of this work included a
`REPORT_TRACK` subworkflow; it duplicated logic the workflow already had and
has been removed rather than left in the tree looking authoritative.

---

## Files

| File | Purpose |
|---|---|
| `modules/local/igv_snapshots.nf` | Breakpoint pages and the clinical SNV page, per sample |
| `modules/local/embed_report_assets.nf` | Inlines each report's dependencies |
| `modules/local/report_zip.nf` | Packages the bundle for release |
| `bin/apply_igv_report_wiring.py` | Applies the workflow and config edits |
| `conf/report.config` | Optional publishDir overrides, not required |

---

## Apply

```bash
cd /goast/mm-awgs-nextflow
python3 bin/apply_igv_report_wiring.py --dry-run
python3 bin/apply_igv_report_wiring.py
```

Five edits: three includes, the IGV stage after `HG38_TRACK`, IGV mixed into
the `ready` gate, the two packaging steps inside the `skip_dashboard` block,
and the parameters in `nextflow.config`. Every anchor is validated before
anything is written; a renamed anchor aborts with nothing touched. Backups per
file, idempotent, `--dry-run` shows the plan.

### The join

```groovy
igv_input = T2T_TRACK.out.mm_annotated_tsv
    .join(T2T_TRACK.out.t2t_bam_bai)
    .join(HG38_TRACK.out.hg38_bam_bai)
    .join(HG38_TRACK.out.v6_report, remainder: true)
    .filter { it[0] != null && it[1] != null }
    .map { meta, mm, tbam, tbai, hbam, hbai, clin ->
        tuple(meta, mm, clin ?: [], tbam, tbai, hbam, hbai)
    }
```

Two decisions worth understanding before this is modified.

**Joined on meta, never combined.** A `combine` would pair every sample's SV
table with every sample's alignments. In a clinical report that means one
patient's breakpoints rendered against another's reads.

**`v6_report` uses `remainder: true`.** A sample with no on-panel clinical
SNVs never emits a clinical table. A plain join would drop that sample from
the IGV stage entirely, taking its translocation pages with it, silently. With
remainder the sample survives and the clinical table arrives as `[]`, which
the process reads as "no somatic snapshots" rather than an error. This is not
hypothetical: one sample in the current cohort has zero on-panel SNVs.

### Publishing

`IGV_SNAPSHOTS` publishes to `${params.outdir}/igv` from a directive in the
module rather than from `conf/modules.config`. That is deliberate:
`REPORT_BUNDLE` scans the published tree, so if the pages are not on disk under
`outdir` before it runs, the report has nothing to link to. Putting the
directive in the module keeps that guarantee with the process that depends on
it. `conf/report.config` carries overrides for anyone who prefers them in
`modules.config`.

### Verify

```bash
nextflow run main.nf -profile docker -stub-run \
  --sample_sheet samplesheets/<sheet>.csv --outdir results_stub

grep -E 'IGV_SNAPSHOTS|REPORT_BUNDLE|DASHBOARD|EMBED_REPORT_ASSETS|REPORT_ZIP' .nextflow.log | tail
```

All five present, `REPORT_ZIP` last. `-resume` is unreliable on this install,
so the stub run is the check rather than a cheap re-run.

---

## Publication

Explicit file list, never `git add -A`.

```bash
cd /goast/mm-awgs-nextflow
git status --short

git add \
  bin/igv_snapshots.py \
  bin/embed_report_assets.py \
  bin/check_igv_pages.py \
  bin/alias_variant_table.py \
  bin/apply_report_rebuild.py \
  bin/apply_igv_report_wiring.py \
  bin/build_report_bundle.sh \
  bin/dashboard_builder/parsers/translocations.py \
  bin/dashboard_builder/parsers/ichor.py \
  bin/dashboard_builder/parsers/qc.py \
  bin/dashboard_builder/templates/translocations_tab.html.j2 \
  bin/dashboard_builder/templates/ichor_tab.html.j2 \
  bin/dashboard_builder/templates/qc_tab.html.j2 \
  bin/dashboard_builder/templates/overview_tab.html.j2 \
  bin/dashboard_builder/build.py \
  bin/dashboard_builder/templates/sample_report.html.j2 \
  tools/make_igv_snapshots.sh \
  bin/make_report_zip.sh \
  modules/local/igv_snapshots.nf \
  modules/local/embed_report_assets.nf \
  modules/local/report_zip.nf \
  conf/report.config \
  workflows/*.nf \
  nextflow.config \
  docs/REPORT_REBUILD.md \
  docs/NEXTFLOW_WIRING.md \
  docs/CHANGES_report_rebuild.md

git status --short
```

```bash
git commit -m "Report rebuild: assay-matched tabs, IGV snapshots, self-contained reports, Nextflow wiring

Replaces report tabs inherited from the hybrid-capture pipeline with ones
matching this assay, and wires the report chain into the workflow.

Dashboard
- QC tab reads QC_ONTARGET products; Picard HsMetrics and fastp removed
- Translocations tab: sortable table, rearrangements by default, paired
  breakpoint IGV viewer
- Copy number tab: ichorCNA figure rasterised for inline display
- Overview reports assay metrics in place of capture metrics
- FLT3 tab and CNV call table removed

IGV
- Per-breakpoint pages plus a manifest, selected by sv_type; rendering the
  whole merged callset produced 3,797 pages and 2.6 GB across three samples
- Clinical SNV page published under the filename the builder resolves

Reports
- Stylesheets, scripts and figures inlined; reports no longer depend on
  their surrounding directory
- Variant tables carry alias columns so the browser resolves every field

Nextflow
- IGV_SNAPSHOTS, EMBED_REPORT_ASSETS, REPORT_ZIP
- IGV joined on meta with remainder on the clinical table, so a sample with
  no on-panel SNVs keeps its translocation pages
- IGV mixed into the existing ordering gate

No variant, gene pair, breakpoint coordinate, FISH finding or expected
karyotype is encoded in any of these files."

git push origin main
```

### Before pushing

```bash
git ls-files | grep -cE '\.bak|nohup|\.preembed'    # want 0
git diff --cached --stat | tail -5
```

Results trees, IGV output, report bundles and archives stay untracked. Confirm
`.gitignore` covers `report_*/`, `results_*/` and `*.zip`.


---

## Scripts must live in bin/

Nextflow adds `<projectDir>/bin` to PATH for every task. `tools/` is never
staged, so a process calling a script from there fails with `command not
found`. `make_report_zip.sh` was moved to `bin/` for this reason; invoke it as
`bin/make_report_zip.sh` when running by hand.

`tools/make_igv_snapshots.sh` stays where it is. It is a standalone wrapper for
results directories that predate the pipeline stage, and is never called from a
process.

### The failure that hid behind this

`build_report_bundle.sh` located `alias_variant_table.py` through `SCRIPT_DIR`.
Run by hand the two sit together; run as a process, the script is staged into a
task directory alone and the sibling is not there. The alias step was guarded
by a `cp` fallback, so `REPORT_BUNDLE` reported success while writing variant
tables with no alias columns, and the Variants tabs came back empty with
nothing in any log to say why.

The helper is now resolved through `SCRIPT_DIR`, then PATH, and the run aborts
if neither finds it. The fallback is gone: a bundle that looks complete and
produces an empty report is worse than a run that stops.

Worth applying the same test to any other script that reaches for a sibling by
path:

```bash
grep -rn 'SCRIPT_DIR\|dirname "\${BASH_SOURCE' bin/*.sh
```
