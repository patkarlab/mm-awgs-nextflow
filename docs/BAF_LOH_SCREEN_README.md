# Panel BAF / LOH screen

`bin/baf_loh_screen.py` — allelic imbalance screen across panel windows, with an
optional copy number join for distinguishing copy-neutral LOH from deletion.

`bin/plot_baf_cn.py` — companion BAF and copy number figures.

---

## What it does

For each window in the panel BED, independently, the screen reads heterozygous
SNV sites from a phased Clair3 VCF and asks whether the B-allele frequency
distribution in that window is consistent with balanced heterozygosity.

Windows are analysed **in isolation**. There is no segmentation, no chaining of
signal between windows, and no assumption that coverage exists between them.
This is deliberate and is the reason the screen works on adaptive sampling data:
coverage forms dense islands over panel targets with very little in between, and
any method that segments across those gaps averages real signal in the islands
against noise in the deserts.

Copy number is not inferred from the BAF. Allelic imbalance is produced by
deletion, by gain, and by copy-neutral LOH alike; separating them requires a
copy number call at the same locus, which is joined in from ichorCNA when
available.

---

## Inputs

| input | notes |
|---|---|
| Panel BED | four columns, chr-named, matching the VCF reference build |
| Phased Clair3 VCF | per-chromosome directory or a single merged file |
| ichorCNA results (optional) | `<dir>/<sample>/ichorcna_out/<sample>.cna.seg` and `.params.txt` |

The BED must be on the same reference build as the VCF. An hg38 VCF with a
T2T-coordinate BED will not error — contig names match — but every window lands
on the wrong sequence.

---

## Usage

Cohort run, with copy number:

```bash
python3 bin/baf_loh_screen.py \
  --bed assets/aWGS_PCN_v7_hg38.bed \
  --sample-map <run>/baf_loh/samples.tsv \
  --ichor-dir <run>/hg38/calls/ichorcna \
  --out <run>/baf_loh/cohort.baf_screen.tsv
```

`samples.tsv` is two tab-separated columns: sample identifier and VCF path. Run
all samples in a single invocation — the cohort het-density baseline requires at
least three samples and is withheld below that.

Single sample, BAF only:

```bash
python3 bin/baf_loh_screen.py \
  --bed assets/aWGS_PCN_v7_hg38.bed \
  --vcf-dir <path to phased VCF dir> --sample <ID> \
  --out <ID>.baf_screen.tsv
```

Figures:

```bash
python3 bin/plot_baf_cn.py \
  --screen <run>/baf_loh/cohort.baf_screen.tsv \
  --bed assets/aWGS_PCN_v7_hg38.bed \
  --sample-map <run>/baf_loh/plot_inputs.tsv \
  --outdir <run>/baf_loh/figures
```

---

## Output columns

| column | meaning |
|---|---|
| `n_het` | usable heterozygous SNV sites in the window |
| `het_per_mb`, `het_density_ratio` | density, and density relative to the cohort median for that same region |
| `median_dp`, `frac_phased` | depth and phased fraction at those sites |
| `median_baf` | unfolded alternate-allele fraction; near 0.5 when balanced |
| `frac_central` | fraction of sites inside the balanced band; the most directly interpretable measure |
| `band_lo`, `band_hi` | the balanced band actually applied, after any ploidy adjustment |
| `baf_deflection` | mean absolute displacement from 0.5 |
| `depletion_score` | observed central density divided by its depth-matched binomial expectation |
| `bimodality` | outer-lobe mass divided by central mass |
| `flag` | `LOH_LIKELY` / `EQUIVOCAL` / `NO_LOH` / `UNASSESSABLE` |
| `cn`, `cn_event`, `cn_call` | copy number at the window and its interpretation |
| `tumour_fraction`, `cn_note` | ichorCNA purity and any caveat on the copy number call |

`LOH_LIKELY` requires **both** central depletion and bimodality. Either alone
gives `EQUIVOCAL`.

`cn_call` values: `CN_LOH` (copy number 2, so imbalance there would be
copy-neutral), `DELETION`, `GAIN`, `GAP` (no ichorCNA bin overlaps and none
could be inferred), `GAP_INFERRED` (inferred from concordant flanking bins),
`NO_CN` (no copy number data supplied).

---

## Design decisions worth knowing

**B-allele frequency is unfolded.** It is computed as `alt / (ref + alt)` across
the full 0–1 range, never folded to `min(ref, alt) / depth`. Folding collapses
the two modes of an LOH distribution onto each other and removes the signal
entirely; it also displaces the median of a balanced region away from 0.5 in a
depth-dependent way, which manufactures apparent deflection everywhere.

**Bimodality is only evaluated at copy number 2.** Above that the balanced state
is itself two-lobed — a normal triploid region sits at 1/3 and 2/3 — so a
two-lobed distribution carries no information about LOH. Where the sample's
modal copy number is not 2, the test is disabled and central depletion alone
downgrades to `EQUIVOCAL` rather than calling LOH on weaker evidence.

**The balanced band follows the sample's modal copy number.** At copy number 3
it shifts to roughly 0.25–0.41 rather than remaining centred on 0.5.

**Het density is corroborating, never primary.** It is confounded by window
size, mappability and coverage, and the wide immunoglobulin-locus windows carry
intrinsically low density in every sample. Comparing each region only against
the cohort median for that same region cancels the confounder.

**Regions below the site minimum report `UNASSESSABLE`, never `NO_LOH`.**
Absence of evidence is reported as such.

**Copy number inference across gaps is conservative.** Where no ichorCNA bin
overlaps a window, the flanking bins are consulted; the inference is taken only
when both flanks exist, agree, and lie within 1.5 Mb, and is labelled
`GAP_INFERRED` rather than presented as a measurement. Inference is suppressed
entirely when the ichorCNA fit is flagged as unreliable.

**ichorCNA fit quality is checked.** A reported tumour fraction of zero, or a
chosen solution whose log-likelihood is below the best candidate in the
solution table, raises a warning that is propagated to every row for that
sample.

---

## Validation

Three samples, one carrying a FISH-confirmed del(17p) in 65% of cells with a
somatic TP53 missense variant at VAF 0.93.

| sample | n_het | frac_central | bimodality | flag | CN | call |
|---|---|---|---|---|---|---|
| A | 736 | 0.694 | 0.44 | NO_LOH | 3 | — |
| B | 693 | 0.651 | 0.54 | NO_LOH | 2 | — |
| C | 160 | 0.062 | 15.0 | LOH_LIKELY | 1 | DELETION |

Sample C reproduced a het-site collapse and BAF bimodality found independently
by manual analysis, and the copy number join correctly resolved the mechanism as
hemizygous deletion rather than copy-neutral LOH. Concordant with FISH,
ichorCNA, and the somatic variant allele fraction.

Panel-wide flag counts and the separation at the validated window were identical
for minimum site depths of 5, 8, 10 and 12.

---

## Limits

**No copy-neutral LOH positive has been observed in the validation cohort.** The
screen has been shown to detect allelic imbalance and to correctly *exclude*
copy-neutral loss when copy number is available. Its behaviour on a true
copy-neutral event is untested.

**The copy number join is weak on adaptive sampling data.** In the validation
cohort no panel window overlapped an ichorCNA bin directly — every call was a
gap or an inference from flanking bins — and two of three ichorCNA fits raised
quality warnings. Treat `cn_call` as advisory. Where a clinical statement
depends on copy number, prefer FISH or an orthogonal assay.

**Genome-wide allele-specific copy number is not achievable from this data.**
Haplotype-resolved copy number callers require phase blocks chained across
continuous coverage; adaptive sampling provides dense islands separated by near
empty regions, and haplotype-specific depth averages toward the mean as a
result. Uniform coverage — for example a shallow whole-genome lane alongside the
enriched run — is the prerequisite, not a different tool.

**For TP53 specifically, variant allele fraction plus copy number is the
stronger route.** A high-VAF somatic TP53 variant with no copy loss implies
copy-neutral LOH directly and needs no BAF analysis. The screen's role there is
corroboration, and it covers the loci where no somatic variant is called.

---

## Environment

Standard library only for the screen; matplotlib additionally for the plots. No
pandas. Runs as-is in the `awgs_sv` conda environment.
