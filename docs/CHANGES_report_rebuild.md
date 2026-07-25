# Report rebuild — change log

The per-case report was a port of the lab's hybrid-capture dashboard. Its tabs
described a different assay: Picard HsMetrics, fastp, per-exon coverage,
FLT3-ITD, CNVkit segment tables. On adaptive-sampling data those render as
dashes or empty panels. This is the record of what changed and why.

## Dashboard

| Tab | Before | After |
|---|---|---|
| Overview | Picard %target-100x, fold-80, AT/GC dropout, low-coverage exons, FLT3 | Panel regions, rearrangement count, tumour fraction, ploidy, clinical variants, reads |
| QC | Picard HsMetrics, fastp, per-exon coverage | mosdepth per-region coverage, read-length and mean-Q plots |
| Variants — Clinical | Loaded 15 rows, showed 3 fields | All fields; GeneBe links resolve |
| Translocations | absent | Sortable table, default sort on supporting reads, paired-breakpoint IGV |
| Copy number | CNVkit segment table and per-gene plots | ichorCNA figure, rasterised, plus fit parameters |
| BAF / LOH | figures broken | figures inlined |
| FLT3 | present | removed |

## Defects found and fixed

**Every tab empty.** The builder passes context as a single `ctx` object; the
new includes used bare names, which Jinja resolves to falsy, so every guard
fell through to its empty state. Includes now resolve under either convention,
and a missing context variable renders a distinct red panel rather than the
same grey box as "parser ran, found nothing".

**Cross-sample lookup.** Two parsers ended with fallbacks that accepted any
matching file in the tree when the sample's own was absent. Inside a
well-formed bundle it never fired, but the failure mode is one sample's report
showing another's data. Both are now bound to the sample identifier with no
fallback.

**IGV rendering the whole callset.** `mm_annotated.tsv` is annotation over the
merged SV callset, not a translocation list. Rendering a page pair per row
produced 3,797 pages and 2.6 GB across three samples. Selection by `sv_type`
brings that to 76 pages for a typical sample. `min_callers` defaults to 1
deliberately: single-caller rearrangements at low support are what this panel
exists to recover.

**Variant browser reading nothing.** The browser resolves columns by exact
name. This pipeline's filter emits lowercase names where the browser expects
capitalised ones; only `REF_COUNT`, `ALT_COUNT` and `Filter` matched, and those
were the only three fields that rendered. The bundle copy now carries alias
columns, and `EXON` is composed from `exon_rank`/`exon_total` since no single
source column holds it.

**Reports depending on their directory.** Stylesheets, scripts and figures were
relative references. They resolved on the server and broke when a report was
copied or extracted elsewhere, silently. All are now inlined.

**Archive defects.** `.preembed` backups shipped inside the clinical archive;
`--light` stripped the IGV tree while leaving its buttons live; a 225 MB
tarball was written alongside a zip of the same tree. Fixed, guarded, and made
opt-in respectively.

## Deliberately not done

**`SomaticSeq_Verdict` and `VariantCaller_Count`** are read by the browser and
drive two filter controls. Both can be populated with `--verdict-from` and
`--caller-count`, and both are off by default: this pipeline runs ClairS-TO,
not SomaticSeq, and a single somatic caller. An empty filter is more honest
than one presenting this pipeline's output under another tool's name or
asserting a caller count that was never measured.

**Mosaic SV track** remains parked. Known clonal Ig translocations sit above
the mosaic VAF floor, so the track adds artefact risk without clinical benefit
in this application.

## Open

- `translocations.tsv` (Ig-aware collapsed, 26 rows) versus `mm_annotated.tsv`
  (annotated full callset, 38 TRA) — the tab renders the second per the
  original specification. Which should be the reporting view is unresolved.
- Clinical and filtered variant tables are byte-identical for the samples
  checked. Either the filter writes the same content to both, or every
  reportable variant is also clinical.
- Wakhan / CN-LOH integration, unchanged from before.


---

# Second pass — fixes found by using the reports

Everything below was found by opening the built reports rather than by testing
the code, which is worth recording: each one produced a report that looked
finished.

## IGV cross-links from variant cards

Every clinical variant card showed "no IGV". The report existed and the hash
router was injected, so the failure was in the lookup that joins them.

`parsers/igv.py` builds its lookup from the igv-reports `tableJson` headers
`CHROM, POSITION, REF, ALT`, keyed as `chr:pos:ref:alt`, and the browser builds
`_igvKey` as `Chr:Start:Ref:Alt`. The sites file written by `igv_snapshots.py`
carried `chrom, start, end` and **no REF or ALT column at all**, so the lookup
came back empty.

The somatic sites file now writes `CHROM, START, END, POSITION, REF, ALT` ahead
of the annotation columns. `POSITION` is the 1-based coordinate; `START` stays
0-based for igv.js. They are deliberately separate columns and should not be
collapsed.

Translocation pages were never affected: the paired viewer resolves through its
own manifest, not this lookup.

## Cohort index

Mean cov, % >= 100x, Fold-80, % dup and FLT3-ITD all came from Picard or a
per-exon coverage table and rendered blank on every row. Replaced with median
depth, regions below 10x, and rearrangement count.

Median rather than mean: the panel mixes focal gene-body windows with megabase
breakpoint windows, so a panel-wide mean is pulled down by the wide ones. The
low-coverage count sits beside it so a good median cannot hide a set of regions
with nothing in them.

## Variant filter thresholds

ALT count offered > 10, > 15, > 20. Against an observed distribution of 63
calls at 1 alt read, 18 at 2 and 8 at >= 3, all three buttons returned nothing.
Now > 1, > 2, > 5.

The Callers filter and its sort entry are removed. The field is
`VariantCaller_Count`, somatic calling here uses one caller, and no option
could ever match. Removed rather than backfilled with a constant: asserting a
caller count that was never measured puts a number in a clinical report that no
part of the pipeline produced. `alias_variant_table.py --caller-count 1`
remains available for anyone who wants the column.

## BAF/LOH tab

`cn_note` repeated the same sentence on all 66 regions, so each row ran to
several lines and the numeric columns became unreadable. The note is now
collapsed to one line and expandable, and when it is identical across every row
it is lifted into a single banner above the table, since in that case it
describes the ichorCNA fit rather than the regions.

A flag filter is built from the values present, with counts, anchored so
`LOH_LIKELY` does not also match a longer value.

Figures gain an include-in-report control, feature-detected against
`window.tspipeReporting`. If that API is absent or named differently the
control does not appear, rather than appearing and silently discarding
selections.

All three are applied to the rendered table by
`templates/baf_loh_enhance.html.j2`, not by rewriting `baf_loh_tab.html.j2`.
The enhancement cannot break the pane and survives changes to it.

## Script resolution under Nextflow

`REPORT_ZIP` called `make_report_zip.sh` from `tools/`, which Nextflow never
stages; only `bin/` is added to PATH. Moved.

`build_report_bundle.sh` located `alias_variant_table.py` through
`SCRIPT_DIR`. Run by hand the two sit together; run as a process the script is
staged into a task directory alone. The alias step had a `cp` fallback, so
`REPORT_BUNDLE` reported success while writing variant tables with no alias
columns and the Variants tabs rendered empty, with nothing in any log to say
why. The helper is now resolved via `SCRIPT_DIR` then PATH and the run aborts
if neither finds it. The fallback is gone.

Any other script reaching for a sibling by path has the same exposure:

```bash
grep -rn 'SCRIPT_DIR' bin/*.sh
```

## Cohort BAF/LOH figures

The pipeline module writes `baf_cn_figures/` where the standalone script writes
`figures/`. The bundle matched only `figures/`, so pipeline runs silently
carried no cohort plots. The directory is now located by content.

## Not fixed, and not a report problem

On-target depth for the 20260720 batch runs 2-6x against a v7 target of
15-20x. At that depth a single alt read is a 30-50% VAF, which is why the
sample with the lowest depth carries the most clinical variants and the highest
mean VAF. No filter threshold recovers information that was not sequenced. This
is a sequencing yield question, not a reporting one.
