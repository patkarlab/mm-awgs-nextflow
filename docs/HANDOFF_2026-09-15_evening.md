# MM aWGS pipeline — handoff, 15 September 2026 (evening revision)

Supersedes the morning version of this file. Read `/mnt/project/project_brief.md`
for background that predates both.

The morning handoff was wrong in two places worth flagging before anything else:
it described `patch_genebaf_bins.py` as "written and tested but not yet applied"
when that script never existed, and its entire ichorCNA section is superseded —
every near-triploid call in the project before this evening was an artefact of
the purity/ploidy fit.

---

## Where things stand

Items 1–4 of the old open list are done and committed. The report's copy-number
tab now reads ichorkaryo's JSON, mmplot draws ichorCNA output, and ichorCNA
itself has been re-fitted against FISH evidence.

Current run resumes session `b079f29f-e83b-447e-b55d-c8836e10c798`. **Always
pass that session id explicitly** — see Infrastructure below; a bare `-resume`
is what filled the disk.

```bash
cd /goast/mm-awgs-nextflow
setsid nextflow run . -profile gandalf \
    -resume b079f29f-e83b-447e-b55d-c8836e10c798 \
    --sample_sheet samplesheets/samplesheet_v7.csv \
    --outdir results_v7 -bg < /dev/null > logs/v7.log 2>&1 &
```

---

## The ichorCNA ploidy finding

This was the substantive result of the day and it invalidates earlier
copy-number output.

Purity and ploidy are degenerate against read depth: a tumour at 67 percent
cellularity with a triploid genome and one at 100 percent with a diploid genome
produce nearly the same log2 profile. ichorCNA was built for cfDNA and its
default normal-fraction grid reflects that — `c(0.5,0.6,0.7,0.8,0.9,0.95)`,
every starting point assuming the sample is at most half tumour. These are
sorted plasma cells at 92–99 percent purity, so the correct answer is a normal
fraction near zero and it was never explored.

With both parameters free, **7 of 15 v7 samples were called near-triploid**. All
seven have cytogenetics, and none supports it:

| sample | FISH centromeres 3, 7, 5, 9, 15 | FISH verdict |
|---|---|---|
| 11F202611457 | x3–4, x3–4, x2, x3–4, x3 | hyperdiploid, trisomies in 20–42% |
| 11F202615394 | all x2, purified CD138+ | diploid, t(11;14) in 78% |
| 11F202617209 | all x2 | normal, 200/200 cells |
| 11F202617530 | all x2, purified CD138+ | diploid, monosomy 13 + 1q gain 95% |
| 11F20265231 | all x2 | diploid, monosomy 13 and 14, 1q amp |
| 11F20265533 | x3, x3, x3, x3, x3 | hyperdiploid, trisomies in 50–65% |
| 12F2023670 | all x2 | "no evidence of hyperdiploidy" |

Five are plainly diploid. The two with trisomies are **hyperdiploid** — a
diploid modal copy number carrying six to eight trisomies — which is not a
near-triploid genome where most of the sequence sits at three. Three reports
were on purified CD138+ plasma cells, removing the usual objection that
unsorted marrow dilutes the signal.

### What was tried, in order

1. **Fix the normal fraction from flow purity.** `--normal c(1 - meta.purity)`
   floored at 0.01, `--estimateNormal FALSE`. This worked — `cna_burden` now
   tracks the flow value on every sample where it previously read 0.668 against
   a flow purity of 0.9968. It did **not** fix ploidy.

2. **Force ploidy.** At pinned cellularity ichorCNA still preferred 3, by 73 log
   units on 11F20265231 (loglik 1373 vs 1300). It wins by widening the states
   that fit badly rather than by moving ploidy: fitted t-precisions on that
   sample were 6.1 at the CN 1 state against 950 at CN 3. A state 150 times
   wider than its neighbour is absorbing those segments, not fitting them. No
   depth-only criterion resolves this.
   `--ploidy c(2)` with `--estimatePloidy FALSE`.

After both, all 15 samples read `baseline_copy_number` 2.

### This is an override, not a correction

ichorCNA was doing what it was asked and the ploidy-3 solution genuinely has the
higher likelihood. The justification is external evidence. With cellularity and
ploidy both pinned there is no free parameter left — copy number is now a direct
reading of log2 against a diploid genome, and ichorCNA is segmenting and
assigning states rather than doing model selection.

`params.ichorcna_estimate_ploidy` exists so a genuinely near-triploid or
near-tetraploid sample can be handled. **That decision needs FISH, not the fit.**

---

## Results after the fix

| sample | ploidy_class | canonical trisomies |
|---|---|---|
| 11F202523808 | non_hyperdiploid | — |
| 11F202610129 | non_hyperdiploid | 19 |
| 11F202611457 | non_hyperdiploid | — |
| 11F202611829 | hyperdiploid | 3, 9 |
| 11F202615394 | non_hyperdiploid | — |
| 11F202617209 | non_hyperdiploid | — |
| 11F202617530 | non_hyperdiploid | — |
| 11F202618760 | non_hyperdiploid | 9 |
| 11F202618931 | non_hyperdiploid | — |
| 11F20265231 | non_hyperdiploid | — |
| 11F20265533 | hyperdiploid | 5, 7, 9, 11, 19 |
| 11F20269040 | non_hyperdiploid | — |
| 11F20269774 | non_hyperdiploid | — |
| 12F2023670 | hyperdiploid | 7, 11 |
| 12F20264455 | hyperdiploid | 3, 7 |

Agreement with FISH where it exists: 265533 hyperdiploid (FISH trisomy 3, 5, 7,
9, 15 in 50–65%), 617209 non_hyperdiploid with nothing (FISH normal on 200
cells), 615394 / 617530 / 265231 all non_hyperdiploid against diploid FISH.

**Two disagreements remain**, both explained and neither a ploidy problem:

`11F202611457` — FISH hyperdiploid, pipeline non_hyperdiploid. The trisomies are
in the data. chr7 segments into 104 bins at log2 0.160 (called CN 2) and 54 bins
at 0.256 (called CN 3); chr9 sits at 0.148 across the whole chromosome. Those
values are exactly what a trisomy at 20–42% clonality predicts. They are lost
twice: ichorCNA discretising 0.160 to CN 2, then `chromosome_copy_numbers`
taking a length-weighted majority that discards the CN 3 portion (chr7 reports
CN 2 at fraction 0.653).

`12F2023670` — pipeline calls hyperdiploid on trisomies 7 and 11; FISH says
none. This is the suspected two-person mixture (open item below). A mixture
produces exactly this: intermediate copy states rounding to 3. **Run
`bin/mmfingerprint.py` before trusting any number from this sample.**

---

## Ploidy is read from the plot

Decision taken this evening: the genome-wide figure makes ploidy and
chromosome-level gains/losses obvious to a reader, and that judgement is more
reliable than the fit. Consequences:

- Do **not** tune `hyperdiploidy()` thresholds until they agree with FISH on
  fifteen samples. That is overfitting against two disagreements, one of which
  may be a contaminated sample.
- `ploidy_class` and the trisomy list should be presented in the report as
  computed-and-to-be-confirmed, not authoritative.
- Report per-chromosome mean log2 alongside the integer call, so a reader sees
  0.160 for what it is.
- The pipeline's job is a figure faithful enough to read. That now requires the
  baseline to be right — the first near-triploid figure looked internally
  coherent and nothing on the page revealed the error.

---

## Majority-vote defect (general, affects three places)

Arms and chromosomes are summarised by length-weighted majority, which discards
partial events by construction.

- **Arms.** 1p on 11F20265231 reads CN 3 while a 35 Mb deletion inside it sits
  at CN 1. `arm_fraction` records the coverage and the report now shows it as an
  "Arm covered" column, amber below 80%.
- **Whole chromosomes.** The 611457 trisomy case above.
- **genebaf arm LOH.** Pools only the panel windows on an arm. On 1p that is
  NRAS and TENT5C — 3 Mb apart and both inside the deleted 35 Mb — so the
  verdict is labelled 1p but describes one segment of it. `--min-arm-windows 2`
  does not catch clustering. **Outstanding fix: record the span the pooled
  windows actually cover alongside each arm verdict.** With BAF scoped to eight
  wide windows there will be more such arms, not fewer.

---

## QC metrics

Detection limits are sound and comparable to mmkaryo's 3.2–15.1%, slightly worse
as expected at 1 Mb bins: **6.5% to 25.7% for a 5 Mb event, median ~9%**.

`sigma_log2_per_bin` ranges 0.036–0.148. 611457 at 0.148 is the outlier and it
is real sub-chromosomal structure at subclonal amplitude, not a bad library —
most of its chromosomes sit at sd 0.055–0.086 and five carry the variance
(chr8 0.497, chr16 0.443, chr18 0.335, chr7 0.325, chr12 0.302). Segmentation on
that sample is normal: a handful of chromosomes split in two or three, same as
265231.

**Overdispersion is measured against the wrong null and should be relabelled or
dropped.** It reads 3.4 to 70 across the cohort, median ~6.6, never near 1. The
metric compares per-bin log2 scatter against Poisson counting noise, but
readCounter counts read starts, ONT read lengths vary widely, and adaptive
sampling rejects non-uniformly — a bin's coverage is not a Poisson count. As
written the card implies something is wrong on every sample. Either relabel it
as excess-over-counting-noise with no implied target, or remove it.

---

## Infrastructure incidents

**`/goast` hit 100%** (7.0 TB, 856 K free). Every MERGE_MINKNOW task exited 1;
the "no .bam files found in barcode04" error for `11F20269774` was the script
failing for want of space, not a bad sample. That sample now completes normally.

**Root cause: bare `-resume` resumes the *last* run.** After a failure that is
the wrong session, so Nextflow finds no usable cache and recomputes from
scratch. `work/` reached 2.3 TB holding roughly 180 realignment BAMs where 30
would do. **Always pass the session id.**

**Three Nextflow cache DBs are corrupt** — `f695da77`, `d24525fc`, `8f4e1cf4`
(the cohort21 and cohort24 sessions). `nextflow clean` cannot read them and will
keep failing on those sessions.

**1.2 TB deleted by timestamp** rather than by session, since clean could not
reach it:

```bash
find work -maxdepth 2 -mindepth 2 -type d ! -newermt "2026-09-15 07:00" -exec rm -rf {} +
```

This was verified safe first — 90 real T2T BAMs and 90 hg38 BAMs were in the
*recent* set against 21 in the old one. It did cost the MERGE_MINKNOW cache,
which `golden_keller` had inherited from an older session, so ~500 GB of merged
BAMs were rewritten.

**Two Nextflow runs were briefly writing to the same `results_v7`** after a
`pkill` that did not take. Resolved; copy-number outputs are unaffected since
both had the patched module.

**Still on disk with no work directories behind them**: `results_cohort21`
(509 GB) and `results_cohort24` (567 GB). Over a terabyte that cannot be
regenerated by resume any more. Decide whether those cohorts are live.

---

## Decisions taken (no code change needed)

**Panel frozen at v7. No v8, no wet-lab changes.** Old item 11 is closed. A
capture BED cannot be applied retrospectively — the molecules in those flanks
were physically ejected, and the only free widening is read spillover of roughly
one fragment length (~10 kb), not 250 kb. The 29.97 Mb figure in the old handoff
was wrong; the correct number for that recipe is **33.86 Mb**, +13.4% target,
which at fixed enrichment costs ~12% depth everywhere.

**v7 is not uniformly flanked.** Only CCND1, CCND2, CCND3, MAFA, MAFB and
TP53+TNFSF12 sit at ±500 kb (1.00 Mb each). MYC is 5.00 Mb, WWOX/MAF 2.45 Mb,
IGH/IGK/IGL 1.31–2.32 Mb. The remaining 57 of 68 regions are gene body ±~50 kb,
median span 184 kb.

**BAF scope narrowed — decided, not yet implemented.** Depth/log2R on all 68
regions; BAF only at the eight wide windows (CCND1/2/3, MAFA, MAFB,
TP53+TNFSF12, MYC, WWOX/MAF). IGH/IGK/IGL excluded as an explicit hard-coded
rule — V(D)J recombination and somatic hypermutation make allele fractions there
clonal, not germline. Per-window allele panels suppressed below usable het
counts; the gate is observed hets, not window width (DEFB1 gives 149 hets in
108 kb, NRAS 11 in 112 kb). **`genebaf.py` still computes allelic state on every
window it is given; a `--baf-regions` argument is outstanding.**

**Open item 8 closed.** The 11F20265231 1p "three-way contradiction" was an
arm-level summary compared against gene-level depth. The CN 1 segment is
chr1:85–120 Mb; NRAS sits at 114.7 Mb and TENT5C at 117.6 Mb, both inside it.
ichorCNA at log2 −0.97, on-target depth at −1.05 and −0.90, and genebaf LOH at
ratio 0.019 all agree.

---

## What was built today

| commit | change |
|---|---|
| `b384974` | genebaf: per-bin depth track at 1 kb; declare bins output |
| `57fe875` | ichorkaryo: make executable |
| `78eebf1` | ichorkaryo: ploidy class, cytoband table, qc block and warnings |
| `60ea070` | Wire cna_seg and readCounter wig through to ichorkaryo |
| `6f0b07e` | mmplot: draw ichorCNA output; baseline-relative lesion classes |
| — | Publish ichorkaryo and genebaf outputs on the hg38 track |
| — | dashboard: read copy number from ichorkaryo JSON |
| — | ichorcna: fix cellularity from flow purity, force diploid ploidy |

**genebaf per-bin depth.** `<prefix>.genebaf.bins.tsv`, ~29,894 rows per sample
at 1 kb, ~1.2 MB. Free: `window_depth()` already called `count_coverage()` and
threw the profile away. Measured 8m10s with and without bins — the cost was
always `count_coverage` over 26 Mb, not the binning. Bins are length-weighted;
the last bin of a window is short, and averaging bins unweighted gives the wrong
window mean.

**ichorkaryo extensions.** `ploidy_class`, `hyperdiploid`, `trisomies`,
`canonical_trisomies`, `estimated_chromosome_count`, `chromosome_copy_number`,
`cytobands`, `qc`, `warnings`. Trisomies are referenced to a constitutional two,
not to the baseline — otherwise a genome whose modal copy number is three
reports no trisomies at all. `detection_limit_5mb` comes from `.cna.seg` per-bin
scatter; `median_reads_per_bin` needs the readCounter WIG, now copied into
`ichorcna_out/`.

**mmplot on ichorCNA.** New `bin/ichorcna_to_mmplot.py`. Three mmplot changes:
copy number read from the `copy_number` column rather than converted from log2
(the conversion only reproduces ichorCNA's numbers when handed ichorCNA's own
fitted parameters, and at flow purity puts CN 1 segments at 1.5); centromeres
derived from the `--cytobands` acen entries rather than a hardcoded T2T table
(chr9 was off by 36 Mb of span on hg38 data); BAF panel dropped when there is no
allele data. Lesion classes are relative to the derived modal baseline, with
sex chromosomes on a constitutional expectation via `--sex`.

**Dashboard.** Parser rewritten against the ichorkaryo JSON; arms and cytobands
come from the JSON rather than sidecar TSVs ichorkaryo never wrote. Template
badge follows `ploidy_class` instead of a boolean. Arm table's blank
detection-limit column replaced by "Arm covered". JSON nulls were rendering as
the literal string "None" — now blank.

**`build_report_bundle.sh`** now also collects `*.ichorkaryo.json`. It still
collects mmkaryo's `karyotype.json`, `cytobands.tsv`, `arms.tsv` and
`segments.baf.tsv`; those go with item 5.

---

## Open items

1. **`bin/mmfingerprint.py` on `12F2023670`** — staged, never run. Sex mismatch,
   suspected two-person mixture, and its hyperdiploid call contradicts FISH. If
   it is a mixture the sample leaves the cohort, which changes the denominator
   of every cohort statistic.
2. **Item 5: retire mmkaryo, mmbaf, COPY_NUMBER subworkflow.** The report
   currently shows ichorkaryo's tables above mmplot figures drawn by mmkaryo
   against T2T — two callers, two references, one tab. Live inconsistency in the
   thing you would hand to someone. Needs `mmplot.nf` rewired onto the converter
   output, the bundle's mmkaryo collectors removed, and the `copy_all` PNG glob
   scoped to a track (`*copy_number*` will match both once mmplot moves).
3. Record the covered span alongside each genebaf arm verdict (see
   Majority-vote defect).
4. Implement `--baf-regions` in `genebaf.py` for the eight-window scope.
5. Relabel or drop the overdispersion QC card.
6. Wire nanomonsv as a confirmation column.
7. Collapse SV mates on MATEID.
8. v4/v5/v6 cohort, separately — the criterion is depth (6.9–10.4x), not panel
   version.
9. **Seven of fifteen samples have no clinical or filtered variant table in the
   bundle** despite CLAIRS_TO and FILTER_V6_REPORT both at 15/15. The bundle
   cannot find `*somaticseq*` for them — likely a filename difference between
   the two variant paths.
10. **Cohort `baf_loh/` missing from the bundle** despite BAF_LOH_SCREEN and
    BAF_CN_PLOTS completing.
11. Untracked and stale: `bin/apply_copy_number_{config,dashboard,wiring}.py`,
    `assets/mm_excluded_junctions.proposed.tsv`.

---

## Things that stay true from the morning handoff

Cohort split by depth; ichorCNA's `genomeBuild`/`genomeStyle` and seqinfo fixes;
`errorStrategy 'ignore'` on ICHORCNA; ichorCNA's "tumour fraction" being
copy-number deviation rather than purity, carried as `cna_burden`; genome-wide
off-target BAF retired; LOH per arm never per gene; nanomonsv's specificity and
its miss on the 2-read t(11;14); SAVANA's `savana_snp_resource = null`; gandalf
resources at 180 cpus / 1200 GB; Ig V-region junction filtering; samplesheets
gitignored; the environment table; the key paths.

One correction: the hg38 panel BED has **66 regions, not 68** — IGH and IGK are
absent, which is deliberate given hg38's `chr14_KI270*` alt contigs sequester
IGH-side reads at MAPQ 0. Neither affects genebaf: IGK was the only 2p window so
2p was never assessable at `--min-arm-windows 2`, chr14 keeps three windows on
14q, and both are Ig loci already excluded from BAF.

---

## File transfer

Files from Claude arrive in `~/inbox/from_claude/`. Files going to Claude are
staged in `~/inbox/to_claude/` and uploaded from there; zip several together
with `zip -j`. Do not use scp.

`~/inbox/from_claude/` holds roughly 90 files going back months, several of
which are the only record of changes not in git. Archive before clearing:

```bash
cd ~/inbox && tar czf from_claude_archive_$(date +%Y%m%d).tgz from_claude/
```
