# Output layout

The pipeline publishes into `${params.outdir}` (default `results/`).

```
results/
├── t2t/                                  # T2T-CHM13 track
│   ├── bams/                             # realigned BAM + bai per sample
│   │   ├── <sample>.t2t.bam
│   │   └── <sample>.t2t.bam.bai
│   └── calls/
│       ├── sniffles/<sample>.sniffles.t2t.vcf.gz(+.tbi)
│       ├── cutesv/<sample>.cutesv.t2t.vcf.gz(+.tbi)
│       ├── severus/<sample>/             # native Severus output dir + normalized <sample>.severus.vcf
│       ├── merged/<sample>.merged.vcf.gz(+.tbi)   # SURVIVOR-merged, SUPP_VEC tagged
│       └── mm_annotated/<sample>.mm_annotated.tsv # headline SV output
│
├── hg38/                                 # hg38 track
│   ├── bams/<sample>.hg38.bam(+.bai)
│   ├── calls/
│   │   ├── clairs_to/<sample>/           # snv_<sample>.vcf.gz, indel_<sample>.vcf.gz + tbi
│   │   ├── ichorcna/<sample>/            # .params.txt, .seg.txt, .cna.seg, <sample>/<sample>_genomeWide.pdf
│   │   └── annotated_clair3/<sample>/    # .pass.vcf.gz, .annotated.vcf.gz, .vep_summary.html, .all_annotated.tsv, .somatic_candidates.tsv
│   └── clair3_phased/<sample>/           # merge_output.vcf.gz, phased_merge_output.vcf.gz + tbi
│
└── pipeline_info/                        # Nextflow report, timeline, DAG, trace
```

## What to look at first

**Per sample:**
1. **`t2t/calls/mm_annotated/<sample>.mm_annotated.tsv`** — every SV that touched
   the v6 panel, with gene calls on both breakpoint sides, known-MM-partner flag,
   supporting caller count, support reads. Headline SV output.
2. **`hg38/calls/annotated_clair3/<sample>/<sample>.somatic_candidates.tsv`** —
   filtered to gnomAD pop AF < 1% and non-synonymous. Driver SNV candidates here.
3. **`hg38/calls/ichorcna/<sample>/`** — read `.params.txt` for tumor fraction +
   ploidy estimate, view `<sample>_genomeWide.pdf` for large-CNV plot. Look for
   1q21 amp, monosomy 13, broad 17p del, hyperdiploidy.

**Per cohort:**
- Per-sample summary aggregation is **not in v0.1** (planned for v0.3).
- Compare ichorCNA `.params.txt` tumor-fraction estimates across samples manually
  for now.

## File-format notes

- **VCFs** are bgzipped (`.vcf.gz`) and tabix-indexed (`.tbi`) throughout.
- **SURVIVOR `merged.vcf.gz`** preserves per-caller support via `SUPP_VEC` INFO
  field. Bit order: `1=Sniffles, 2=CuteSV, 3=Severus`. The MM annotation step
  decodes this into the `callers` column of the TSV.
- **ClairS-TO** produces two VCFs per sample: `snv_<sample>.vcf.gz` and
  `indel_<sample>.vcf.gz`. PoN-filtered NonSomatic variants are tagged in
  the FILTER column.
- **ichorCNA `.cna.seg`** is `chr start end ... copy.number logR.copyNumber subclone.status`.
- **VEP TSV (`*.all_annotated.tsv`)** has 19 columns; HGVS c-dot and p-dot are
  intentionally omitted (multi-allelic crash workaround per step 12b).
- **Clair3 `phased_merge_output.vcf.gz`** is the phased germline VCF. Reserved
  as input for Wakhan in v0.3.
