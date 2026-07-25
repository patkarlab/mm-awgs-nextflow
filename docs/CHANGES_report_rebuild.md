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
