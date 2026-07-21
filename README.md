# mm-awgs-nextflow

**Adaptive whole-genome sequencing (Oxford Nanopore) for multiple myeloma —
a dual-reference, multi-caller pipeline that uses T2T-CHM13 for IGH
translocation detection and hg38 for SNV / indel / CNV / phased germline
calling, ported from a validated bash pipeline to Nextflow.**

Status: `v0.1-dev` — early development. See the [Scope](#v01-scope-explicit)
and [Roadmap](#roadmap) sections for what's in this version and what isn't.

## Why dual reference

| Reference | Used for | Why |
|---|---|---|
| **T2T-CHM13v2.0** | Structural variants (SV calling, IGH translocation detection) | hg38's chr14 alt-contigs (`chr14_KI270*`) sequester IGH-side reads (MAPQ 0 ~95% on hg38 → <1% on T2T). The T2T-CHM13 IGH region is fully assembled and unambiguous. |
| **hg38** | SNV, indel, large CNV, phased germline, VEP annotation | Mature tooling. ClairS-TO, ichorCNA, Clair3, VEP, and Wakhan all ship with hg38-coordinate reference panels (G1000, gnomAD, ClinVar, GC/mappability wigs, PoNs). Re-training all of those for T2T is out of scope. |

Both tracks consume the same input MinKNOW BAM. They run in parallel — neither
depends on the other.

## What this pipeline does (v0.1)

```
sample sheet
    │
    ├── T2T track ────► Sniffles2 + CuteSV + Severus
    │                       │
    │                       └─► SURVIVOR merge ─► MM partner-pair annotation
    │
    └── hg38 track ───┬─► ClairS-TO (somatic SNV/indel)
                      ├─► ichorCNA (large CNV from off-target reads, per-sample)
                      └─► Clair3 phased germline ─► VEP v113 annotation
                                                       │
                                                       └─► somatic candidates TSV
```

Per-sample outputs:
- **T2T**: three per-caller SV VCFs, SURVIVOR-merged VCF with per-caller `SUPP_VEC`, MM-annotated TSV with gene-pair calls and known partner flags
- **hg38**: ClairS-TO somatic SNV+indel VCFs, ichorCNA segment file + tumor fraction + plots, Clair3 phased germline VCF, VEP annotation TSV with somatic candidates (gnomAD AF < 1% + non-synonymous filter)

## Quick start

```bash
git clone https://github.com/<org>/mm-awgs-nextflow
cd mm-awgs-nextflow

# Smoke test (no real data needed)
nextflow run main.nf -profile stub,conda -stub-run

# Real run on gandalf
cp assets/sample_sheet_template.csv my_samples.csv
$EDITOR my_samples.csv   # fill in real sample IDs + MinKNOW BAM paths

nextflow run main.nf \
    -profile conda,docker,gandalf \
    --sample_sheet my_samples.csv \
    --outdir results/
```

The default config points all reference paths at the production gandalf
locations. Override any with `--<name> <value>` — see `--help`.

## v0.1 scope (explicit)

**In v0.1:**
- T2T realignment + three SV callers + SURVIVOR merge + MM annotation
- hg38 realignment + ClairS-TO + ichorCNA (per-sample) + Clair3 phased + VEP
- Single-tree output (no parallel `_18hrs/` tree)
- Conda + docker mixed environment matching the bash production pipeline

**NOT in v0.1** (in the bash production pipeline but deferred):
- SV per-caller read-support enrichment (`enrich_sv_with_support.py`)
- T2T coverage QC (`03_coverage.sh`)
- DeepVariant as second SNV caller
- Clair3 ∩ DeepVariant consensus at relaxed thresholds
- mpileup validation gates (strict + relaxed)
- TP53 focal LOH (custom AF-distribution test, step 13)
- Wakhan haplotype CN/LOH (step 14)
- IGV report generation
- Per-sample / cohort summary TSV assembly

These layer on as v0.2 and v0.3.

## Known v0.1 limitations

1. **The v6 hg38 panel BED does not exist yet.** ClairS-TO uses the v5 hg38
   BED (matches current production bash). Building a v6 hg38 BED by lifting
   over the 10 new v6 driver regions (CrossMap T2T → hg38) is tracked as a
   v0.2 work item.
2. **The placeholder `assets/aWGS_MMfocused_v5_hg38.bed` is not the real BED.**
   Before first real run, copy the real one from gandalf:
   ```
   cp /goast/nikhil_awgs_testing/hg38/beds/aWGS_MMfocused_v5_hg38.bed \
      assets/aWGS_MMfocused_v5_hg38.bed
   ```
3. **Reference paths default to gandalf locations.** Running on another
   server requires overriding `--t2t_fasta`, `--hg38_fasta`, `--severus_pon`,
   `--ichorcna_env_prefix`, `--vep_cache_dir`, `--clair3_model_dir`, etc.
4. **`/home/hemat` on gandalf is 98% full.** The `gandalf` profile sets
   `workDir = /goast/nikhil_awgs_testing/nextflow_work` to keep work staging
   off `/home`.

## Validation history

The bash pipeline this Nextflow port is based on was validated against
orthogonal clinical assays; cohort-level validation details are maintained
outside this repository. Bit-for-bit validation of the Nextflow port against
the bash output is pending.

## Panel design

The v6 T2T panel covers **24.16 Mb across 38 regions (0.776% of T2T-CHM13v2.0)**:
- 28 v5-locked regions for canonical and rare IGH partners, light-chain
  loci, MYC, and high-risk-prognosis genes (TP53, CDKN2C, RB1)
- 10 additional MM driver genes added in v6 (DIS3, TRAF3, PRDM1, ATM,
  CYLD, H1-4, MAX, EGR1, LTB, ATR)
- TP53 trimmed from ±1 Mb to ±500 kb

Full region table, MM relevance, and PMIDs in
[`assets/aWGS_MMfocused_v6_README.md`](assets/aWGS_MMfocused_v6_README.md).
The panel is reproducibly buildable from the v5 panel + the NCBI T2T-CHM13v2.0
RefSeq GFF using `bin/build_v6_panel.py`.

## Caller tuning (validated v3-era defaults)

| Caller     | Key flags                                                       |
|------------|-----------------------------------------------------------------|
| Sniffles2  | `--minsupport 2 --mapq 0 --min-alignment-length 300 --output-rnames` |
| CuteSV     | `--min_support 2 --min_mapq 0 --min_size 30 --min_read_len 300` |
| Severus    | `--min-mapq 10 --min-support 2` + T2T VNTR + T2T PoN            |
| SURVIVOR   | `500 1 1 1 0 30` (500bp tolerance, 1-caller min, type+strand match) |
| ClairS-TO  | `--platform ont_r10_dorado_sup_5khz_ssrs --disable_verdict`     |
| ichorCNA   | bins 1 Mb, ploidy `c(2,3)`, normal grid `c(0.5..0.95)`, panel-as-centromere mask |
| Clair3     | `--platform ont --enable_phasing --longphase_for_phasing`, model `r1041_e82_400bps_hac_v520` |
| VEP        | v113 merged cache, NO `--hgvs` (multi-allelic crash workaround) |

All overridable via `--<flag> <value>`.

## Environment

Each Nextflow process picks the right conda env or docker image (see
`conf/envs.config`):

| Step | Mechanism |
|---|---|
| All T2T-track tools, REALIGN_HG38, helper python | conda `awgs_sv` |
| ClairS-TO | docker `hkubal/clairs-to:latest` |
| ichorCNA | native conda env `ichorCNA` (absolute paths, no activation) |
| Clair3 phased | docker `hkubal/clair3:latest` |
| VEP | docker `ensemblorg/ensembl-vep:release_113.0` |

## Roadmap

- [x] v0.1: T2T SV ensemble + hg38 SNV/CNV/phased germline + VEP annotation
- [ ] v0.2: v6 hg38 BED (CrossMap), DeepVariant, Clair3 ∩ DV consensus, mpileup gates
- [ ] v0.3: Wakhan CN-LOH, IGV reports, sample/cohort summary assembly
- [ ] v1.0: Validation on n ≥ 20 cohort, peer-reviewed publication

## License

MIT. See [LICENSE](LICENSE).
