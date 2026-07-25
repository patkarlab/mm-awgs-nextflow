# Report rebuild — installation and status

Rebuilds the per-case report around the artefacts this pipeline actually
produces, removing the tabs inherited from the hybrid-capture pipeline the
dashboard builder was ported from.

Everything here runs in `awgs_sv`, which was confirmed to hold both
`create_report` and `jinja2`/`pandas`. No second environment is needed.

---

## What is in this archive

| File | Destination | Purpose |
|---|---|---|
| `bin/dashboard_builder/parsers/translocations.py` | same path in repo | Reads `*.mm_annotated.tsv`, types the support columns as numeric, resolves IGV pages per event |
| `bin/dashboard_builder/parsers/ichor.py` | same path in repo | Reads the ichorCNA figure and `params.txt`, rasterises the PDF for inline display |
| `bin/dashboard_builder/parsers/qc.py` | same path in repo | Reads QC_ONTARGET plots and per-region coverage; computes read N50 and mean Q |
| `bin/dashboard_builder/templates/translocations_tab.html.j2` | same path in repo | Sortable SV table plus the paired-breakpoint IGV viewer |
| `bin/dashboard_builder/templates/ichor_tab.html.j2` | same path in repo | Copy-number tab, figure only, no call table |
| `bin/dashboard_builder/templates/qc_tab.html.j2` | same path in repo | Adaptive-sampling QC tab |
| `bin/igv_snapshots.py` | same path in repo | Generates IGV pages for clinical SNVs and for each translocation breakpoint |
| `bin/build_report_bundle.sh` | replaces existing | Now also collects `igv/` and `baf_loh/` |
| `bin/apply_report_rebuild.py` | same path in repo | Idempotent patch wiring the above into `build.py` and `sample_report.html.j2` |
| `tools/make_igv_snapshots.sh` | same path in repo | Standalone cohort wrapper, for runs that predate IGV being wired in |
| `bin/make_report_zip.sh` | same path in repo | Packages a built bundle into a zip in `~/inbox/from_claude/` |

---

## Install

```bash
cd /goast/mm-awgs-nextflow
source /home/hemat/anaconda3/etc/profile.d/conda.sh && conda activate awgs_sv

# 1. Unpack (files land on their repo paths)
unzip -o ~/inbox/from_claude/mm_report_rebuild.zip -d .
chmod +x bin/igv_snapshots.py bin/build_report_bundle.sh tools/make_igv_snapshots.sh

# 2. Inspect the patch before it touches anything
python3 bin/apply_report_rebuild.py --dry-run

# 3. Apply
python3 bin/apply_report_rebuild.py
```

The dry run prints every block it would delete in full. Read the QC pane and
CNV pane blocks before applying: anything in them worth keeping needs to be
moved into the new includes first.

The patch aborts without writing if any anchor is missing, writes a timestamped
`.bak` next to each file it edits, and refuses to apply twice.

---

## Run

```bash
RUN=results_v7_20260713_24h

# 1. IGV pages. Check what it resolves before rendering anything.
DRY_RUN=1 tools/make_igv_snapshots.sh "$RUN"
tools/make_igv_snapshots.sh "$RUN"

# 2. Bundle, now including igv/ and baf_loh/
bin/build_report_bundle.sh "$RUN" report_v7_20260713_24h

# 3. Dashboard
python3 bin/dashboard_builder/build.py report_v7_20260713_24h
```

Step 1 is the slow one: two `create_report` invocations per event, each
embedding a BAM slice. Run it detached if the cohort is large.

```bash
setsid tools/make_igv_snapshots.sh "$RUN" > igv_snapshots.log 2>&1 < /dev/null &
```

---

## How the paired IGV view works

`create_report` renders one locus per page. Rather than post-processing its
embedded session JSON to force igv.js into multi-locus mode — which would tie
the report to an igv-reports version whose feature set is already known to be
mislabelled — `igv_snapshots.py` writes one standalone page per breakpoint and
a manifest linking them:

```
igv/translocations/<event_id>.A.html
igv/translocations/<event_id>.B.html
igv/translocations/<sample>.translocations.manifest.json
```

The Translocations tab loads the two pages into two iframes side by side.
Nothing is shared between the frames, so there is no `srcdoc` escaping, no
cross-frame scripting, and no dependency on internal structure. Each page's
variant table carries a `partner_locus` and `partner_gene` column so a page
still makes sense opened on its own.

Events with no manifest entry render without an IGV button rather than with a
button that leads nowhere.

---

## Verified before delivery

Parsers, templates and the patch were exercised against fixtures built to the
exact schemas from the run directory: the 21-column `mm_annotated.tsv`, the
28-column `clinical.tsv`, a five-column `region_coverage.tsv`, a 50,000-row
`readlen_qscore.tsv`, and a real three-page PDF.

- All three parsers return correct structures; `known_freq` values such as
  `<1%` and `5-10%` sort by magnitude rather than as text
- Templates render and produce balanced HTML
- `igv_snapshots.py` handles both modes, the zero-row case (placeholder page,
  exit 0), and unparseable coordinates
- The patch is idempotent, aborts cleanly on a missing anchor, and produces
  syntactically valid Python and a renderable template

`create_report` itself was stubbed, since igv-reports is not installed in the
build environment. The invocation is echoed in full on failure.

---

## Still open

**Needs `dashboard_builder_src.tar.gz`.** These four depend on files not yet
seen:

1. **Variants — Clinical with GeneBe.** `build.py` already calls
   `p_genebe.annotate(clinical_rows=...)` behind an `annotate_genebe` flag, so
   this is enabling the flag and rendering the returned columns, not new
   plumbing. Needs `parsers/genebe.py` and the clinical tab markup.
2. **IGV tab pointing at the SNV report.** `build.py` discovers an IGV file and
   sets `ctx["files"]["igv_report"]`, then patches its hash router. The somatic
   page must land where that discovery looks. Needs the discovery block and
   `parsers/igv.py`.
3. **Reporting tab.** Untouched by this patch; the include/exclude checkbox
   flow needs checking against the new Translocations table.
4. **Variants — All Filtered.** Left in place. Confirm whether it should stay.

**Not blocking, worth noting.**

`known_mm_pair` and `known_freq` are populating correctly — the cohort check
returned canonical pairs with their reference frequencies. The token-set
intersection fix worked; that item can come off the pending list.

IGV is still not a pipeline stage. `modules/local/igv_report.nf` does not exist
in this repo and `nextflow.config` has no `igv` parameters, so IGV output will
keep needing the standalone wrapper until a module is written and wired into
the subworkflows. That wiring needs the subworkflow files.


---

## Correction: IGV event selection (2026-07-25)

The first version of `igv_snapshots.py` rendered a breakpoint page pair for
every row of `mm_annotated.tsv`. That table is annotation over the whole merged
SV callset, not a translocation list. On the 20260713 run it produced 3,797
pages and 2.6 GB of embedded alignment data across three samples, because the
row counts break down as:

| sv_type | rows (12F20264455) |
|---|---|
| DEL | 421 |
| INS | 220 |
| TRA | 38 |
| INV | 6 |
| DUP | 6 |

`igv_snapshots.py` now selects events before rendering. Defaults: `sv_type`
in `TRA`, interchromosomal only, `--min-callers 1`, `--max-events 200`. On the
same sample that is 38 events and 76 pages. The selection is printed before
rendering starts, so the cost is visible up front.

`--min-callers` defaults to 1 on purpose. Single-caller rearrangements at low
read support are the calls this panel exists to recover; the FISH-concordant
t(11;14) results in this project came from CuteSV alone at RE=2. Raising the
floor is a per-run decision, not a default.

The Translocations tab renders every row of the annotated table but filters to
rearrangements by default, with a switch to reveal the other SV classes. Rows
without breakpoint pages render without an IGV button rather than with one that
leads nowhere.

### Regenerating after the fix

```bash
rm -rf results_v7_20260713_24h/igv
tools/make_igv_snapshots.sh results_v7_20260713_24h
bin/build_report_bundle.sh results_v7_20260713_24h report_v7_20260713_24h
python3 bin/dashboard_builder/build.py report_v7_20260713_24h
bin/make_report_zip.sh report_v7_20260713_24h --light
```

### Open question: which table the tab should render

`mm_annotated.tsv` holds 691 rows including 38 TRA. `translocations.tsv` holds
26 rows for the same sample. The second is the Ig-aware collapsed set, so it is
the more clinically meaningful list, but the two are not the same events and
the identifiers may not correspond. Worth resolving before the tab is treated
as the reporting view.


---

## Fix: empty tabs (2026-07-25)

Every new tab rendered its empty state on the first real build. The cause was a
context-naming mismatch, not missing data.

`build.py` passes the per-sample context as a single `ctx` object; the existing
panes address everything through it (`ctx.flt3`, `ctx.cnv.clinical_table`,
`ctx.files.fastp`). The new includes used bare `qc`, `ichor` and
`translocations`, which are undefined under that convention. Jinja resolves an
undefined name to a falsy value, so `{% if not qc %}` was always true and every
guard fell through to the empty branch.

Three changes:

1. **Each include resolves its context under either convention**, checking the
   bare name first and falling back to `ctx.<name>`. The includes now work
   whichever way the builder renders them.
2. **Parsers always return a dict**, carrying `found`, `reason` and `searched`.
   A `None` return was indistinguishable from an unbound template variable,
   which is precisely what made this hard to read off a rendered page. The
   empty state now prints the directory that was searched.
3. **A missing context variable renders a distinct red panel**, not the same
   grey empty state as "parser ran, found nothing". Those are different
   failures and should never look alike.

### Cross-sample lookup, found during the same fix

`_find_table` ended with a fallback that accepted any `*.mm_annotated.tsv` in
the tree when the sample's own file was absent, and `_qc_dir` accepted any
directory named `qc`. Inside a well-formed bundle each sample directory holds
only its own files, so this never fired in practice, but the failure mode is
one sample's report displaying another sample's rearrangements. Both lookups
are now bound to the sample identifier with no fallback: a directory that does
not hold this sample's file renders empty.

### Variant tables in the bundle

`build_report_bundle.sh` renamed the SNV tables to `<sample>.clinical.tsv`,
stripping the upstream `.somatic_candidates.v6_clinical.tsv` suffix. If the
builder discovers those tables by glob, renaming them is enough to make
discovery miss, which would leave the Variants tabs empty for the same class of
reason. Both files are now copied twice, once normalised and once under the
original name.

---

## Report zip location

`bin/make_report_zip.sh` writes the archive to the directory the command was
run from. Nothing is hardcoded to any inbox or site path.

```bash
cd /goast/mm-awgs-nextflow
bin/make_report_zip.sh report_v7_20260713_24h --light
# -> /goast/mm-awgs-nextflow/report_v7_20260713_24h_light.zip
```

`--out` overrides the destination when it is needed. A relative `--out` is
resolved against the invocation directory, not against the bundle.

The archive is checked for internal path length before it is handed over.
Windows Explorer still refuses to extract paths beyond 260 characters, and the
bundle nests `<sample>/igv/translocations/<event_id>.A.html`; the script warns
when the margin is thin so the failure happens on the server rather than on a
reporting workstation.


---

## Self-contained reports (2026-07-25)

Relative references were the common cause behind the BAF/LOH figures showing
as broken-image icons and the styling not travelling with a copied report. The
builder writes `assets/css/...`, `baf_loh/figures/....png` and similar as paths
relative to the bundle; they resolve on the server and break as soon as a file
is moved, mailed, or extracted to a different depth.

`bin/embed_report_assets.py` rewrites each generated report so every local
dependency is carried inside it:

- stylesheets become inline `<style>`, including any `url()` they reference
- scripts become inline `<script>`, with `</script>` escaped so a minified
  library cannot terminate its own block
- images become base64 data URIs
- absolute URLs, existing data URIs and page anchors are left alone

IGV breakpoint pages are deliberately **not** inlined. Each already carries its
own alignment slice, so folding them in would multiply the report by the number
of events. They travel beside it, and the script reports how many iframe
targets are present versus missing.

Anything referenced but absent is listed by path at the end of the run. That
list is the bundle step's to-do: a file named there was never collected.

Run it after the dashboard build, before zipping:

```bash
python3 bin/dashboard_builder/build.py report_v7_20260713_24h
python3 bin/embed_report_assets.py report_v7_20260713_24h
bin/make_report_zip.sh report_v7_20260713_24h
```

A `.preembed` copy of each rewritten file is kept; `--no-backup` skips it.

### Other fixes in the same pass

**QC read summary.** `readlen_qscore.tsv` is a one-row summary, not a per-read
table. The parser was consuming that row as a single read, which is why the
card read `Reads 1` with N50 equal to the mean. A one-row file is now detected
and its values displayed directly.

**Overview tab.** Replaced. The inherited panel reported Picard
percent-target-100x, fold-80 penalty, AT/GC dropout, low-coverage exons and
FLT3-ITD, none of which this assay produces, so all of them rendered as dashes.
The replacement reports on-target regions, rearrangement count, tumour fraction
and ploidy, clinical variant count and read count, each read from a parser
result.

**Variants tabs.** The tab text names the file the builder looks for
(`somaticseq_clinical_final.tsv`). Normalising the bundle filenames to
`<sample>.clinical.tsv` is what made discovery miss. The bundle now also writes
aliased copies under the names build.py actually looks for, in both the sample
root and `snv/`:

```
<sample>_somaticseq_clinical_final.tsv
<sample>_somaticseq_filtered.tsv
```

Taken from the builder's own discovery patterns rather than guessed. Override
with `SNV_ALIAS_CLINICAL` / `SNV_ALIAS_FILTERED` if they ever change.

Two other patterns appear in build.py and are deliberately not satisfied:
`<sample>_exon_coverage.tsv` and `<sample>_flt3_consensus.tsv`. Both belong to
the hybrid-capture pipeline. The FLT3 tab is removed, and per-exon coverage is
replaced by per-region coverage from mosdepth, which is the meaningful unit for
a panel of breakpoint windows rather than exons.

**`--light` guard.** `--light` strips `igv/`, but the reports were built
against a bundle that had it, so the IGV buttons stayed live and opened empty
panes, which is what the broken-file icons in the paired viewer were. The zip
tool now refuses `--light` when the reports contain IGV controls, unless
`--force`.


---

## Diagnosing an empty paired-breakpoint viewer

Three different failures put an empty pane on screen and they are
indistinguishable by looking at it:

1. the page is not in the archive (stripped by `--light`, or never collected)
2. the page is there but carries no pileup (no reads at that locus, or a track
   that failed to attach)
3. the page is there and complete, but the report cannot reach it

`bin/check_igv_pages.py` separates them:

```bash
python3 bin/check_igv_pages.py report_v7_20260713_24h
```

It reads each manifest, resolves every referenced page, and reports present /
missing / thin per sample. A working igv-reports page is a few hundred
kilobytes because the read pileup is embedded in it as base64; anything much
smaller opened without error and has nothing to show. The script prints a
one-line diagnosis at the end pointing at which of the three cases applies.

### Removing case 3 entirely

```bash
python3 bin/embed_report_assets.py report_v7_20260713_24h --embed-igv 12
```

This inlines the breakpoint pages for the first N IGV buttons as `data:` URIs,
so the viewer works from the report file alone with no sibling directory to
lose. The table is sorted by supporting reads, so those N are the
best-evidenced events. Buttons past the limit keep their relative paths and
still work while `igv/` is alongside.

Budget roughly the size of 2N igv-reports pages, typically a few hundred
kilobytes each. `--embed-igv 12` on a three-sample cohort adds on the order of
tens of megabytes and makes the top events portable anywhere.


---

## Routine run

```bash
cd /goast/mm-awgs-nextflow
RUN=results_v7_20260713_24h
BUNDLE=report_v7_20260713_24h

tools/make_igv_snapshots.sh "$RUN"
bin/build_report_bundle.sh "$RUN" "$BUNDLE"
python3 bin/dashboard_builder/build.py "$BUNDLE"
python3 bin/embed_report_assets.py "$BUNDLE"
bin/make_report_zip.sh "$BUNDLE"
```

Open `<BUNDLE>/cohort_index.html`. The archive lands beside it as
`<BUNDLE>.zip` in whatever directory the command ran from.

### On --embed-igv

Leave it off for routine work. The full archive already carries `igv/`, so the
pages load from there and each report stays around 2.5 MB, which opens
instantly.

`--embed-igv N` exists for one case: sending a single report HTML detached from
its folder. It inlines the breakpoint pages for the top N events, so those work
standalone. The cost is real -- `--embed-igv 12` took reports from 2.6 MB to
22 MB, and a 22 MB HTML is noticeably slow to open in Chrome. Use a small N, or
none.

### What the embedder's closing report means

- **Assets referenced but not found** -- files the bundle step did not collect.
  Actionable.
- **References that resolve to a directory** -- a malformed src or href in a
  template. Harmless to rendering, but it belongs in the template.
- **IGV pages referenced but not present** -- the paired viewer will show empty
  panes for those events.

An empty list under all three is the signal that the bundle is complete.


---

## Clinical variants tab: loaded but blank (2026-07-25)

Fifteen rows displayed with almost every field empty. The file was found and
parsed; the columns were not.

The variant browser resolves columns by exact name. This pipeline's filter
emits mostly lowercase names, while the browser was written against a table
using capitalised ones:

| filter output | browser expects |
|---|---|
| `gene` | `Gene` |
| `consequence` | `Consequence` |
| `impact` | `IMPACT` |
| `tumor_af_pct` | a VAF percentage column |
| `REF_COUNT` | `REF_COUNT` |
| `ALT_COUNT` | `ALT_COUNT` |
| `Filter` | `Filter` |

The last three match in both, and those were exactly the three fields that
rendered. Everything above them came back empty, which reads on screen like a
missing file but is a naming mismatch on a file that loaded fine.

`bin/alias_variant_table.py` writes the bundle copy with a capitalised
duplicate appended for each column, originals preserved, so this pipeline's own
parsers and the browser both find what they read. No value is transformed or
invented; an alias column is a verbatim copy of its source, and an alias whose
source is absent is not written.

### Redundant files in snv/

The previous version wrote six files there, all byte-identical: a normalised
name, the upstream name, and an alias, for each of clinical and filtered. That
is now two, one per class, under the names build.py discovers. The aliased file
is a superset of the original, so nothing is lost by dropping the duplicates;
the untouched originals remain in the results tree.

Worth checking upstream: for 11F20265231 the clinical and filtered tables are
byte-identical. Either the filter writes the same content to both outputs, or
every reportable variant in this sample is also clinical. The second is
plausible at 15 variants, but if the two files are always identical the filter
is not doing what its two output names imply.

### IGV tab

The tab reads a single page discovered by filename, and that name is not
knowable from outside build.py. The somatic page is now also copied to the
sample root under two candidate names. Pin the right one with `IGV_ALIAS` once
it is known:

```bash
grep -oE "'[^']*\.html'" bin/dashboard_builder/build.py | sort -u
IGV_ALIAS=<name> bin/build_report_bundle.sh <results> <bundle>
```
