# Development

## Repository layout

```
mm-awgs-nextflow/
├── main.nf                    # entry point, help, validation
├── nextflow.config            # params, profiles, manifest
├── conf/
│   ├── base.config            # per-label resource defaults
│   ├── envs.config            # per-process conda/docker assignments
│   ├── modules.config         # per-process publishDir paths
│   ├── gandalf.config         # gandalf-specific local executor
│   └── test.config            # tiny smoke-test profile
├── workflows/mm_awgs.nf       # top-level workflow
├── subworkflows/local/
│   ├── prepare_input.nf
│   ├── t2t_track.nf           # T2T SV pipeline
│   └── hg38_track.nf          # hg38 SNV/CNV/germline pipeline
├── modules/local/             # one process per file
│   ├── realign_t2t.nf, realign_hg38.nf
│   ├── samtools_index.nf      # SAMTOOLS_INDEX_T2T + SAMTOOLS_INDEX_HG38
│   ├── sniffles.nf, cutesv.nf, severus.nf
│   ├── survivor_merge.nf, annotate_mm_translocations.nf
│   ├── clairs_to.nf, ichorcna.nf
│   └── clair3_phased.nf, vep_annotate_clair3.nf
├── bin/                       # helpers staged onto PATH inside processes
│   ├── annotate_mm_translocations.py
│   ├── build_v6_panel.py
│   └── qc_v6_cohort.sh
├── assets/                    # static inputs committed to git
│   ├── aWGS_MMfocused_v6_t2t_chr.bed
│   ├── aWGS_MMfocused_v6_t2t_NC.bed
│   ├── aWGS_MMfocused_v6_README.md
│   ├── aWGS_MMfocused_v5_hg38.bed       (PLACEHOLDER; see assets/aWGS_MMfocused_v5_hg38.bed)
│   ├── mm_translocation_dictionary.tsv
│   └── sample_sheet_template.csv
├── tests/                     # synthetic inputs for stub-run / CI
├── docs/                      # this directory
└── .github/                   # CI workflows + issue templates
```

Modules are intentionally thin. Each `.nf` declares exactly one `process`
with a `script:` block (real tool invocation) and a `stub:` block (touch
expected outputs). Adding a new tool = a new module file + an include in
the appropriate subworkflow.

## How a Nextflow port differs from the bash

The production bash scripts (under `/goast/nikhil_awgs_testing/scripts/v6_hg38/`
and `/goast/nikhil_awgs_testing/t2t/scripts/`) each `source paths.sh` and
`conda activate awgs_sv` at the top. Nextflow modules don't `source` external
files. Instead:

- All path constants live in `params.*` in `nextflow.config`.
- Env activation is declarative via the `conda` or `container` process directives
  in `conf/envs.config`.
- A few processes use absolute paths to env binaries instead of activation
  (`ICHORCNA` uses `${params.ichorcna_env_prefix}/bin/readCounter` etc.) —
  this matches the production bash's `11_ichorcna_hg38.sh` choice.

## Running the stub workflow

```bash
nextflow run main.nf -profile stub,conda -stub-run
```

Walks the whole DAG using each module's `stub:` block. Fast end-to-end check
that channel wiring works without needing real BAMs.

## Adding a new tool

1. Create `modules/local/<tool>.nf` — copy the structure of `modules/local/sniffles.nf`
2. Include it from the relevant subworkflow under `subworkflows/local/`
3. Add any new parameters to `nextflow.config` under `params { ... }`
4. Add publishDir to `conf/modules.config`
5. Add conda/container directive to `conf/envs.config`
6. Update `docs/usage.md` for new flags and `docs/output.md` for new outputs
7. Run `-profile stub -stub-run` to verify wiring

## Validating against the bash production output

Before merging changes to a module that already has a production bash equivalent:

1. Run the bash production script on one sample. Note all output files + sizes.
2. Run the Nextflow port on the same sample.
3. Compare outputs file-by-file. Acceptable differences:
   - Timestamps inside VCF headers
   - Reordered records within position-sort ties
   - Temp-file paths in log output
4. NOT acceptable:
   - Different variant counts
   - Different PASS/FAIL filter assignments
   - Different INFO field values (SUPP_VEC, RE, etc.)
   - Missing output files

## Sample-ID and PHI policy

- No patient names, MRNs, FISH report contents, or other PHI in any committed file.
- Sample IDs in `tests/` and `assets/sample_sheet_template.csv` are synthetic
  (`SAMPLE_01`, etc.).
- Real sequencing IDs only flow through the pipeline at runtime via the user's
  sample sheet, which is not committed.
- The `.gitignore` blocks `*.bam`, `*.vcf*`, `patient_*`, `fish_*`, etc. by
  default. Do not weaken these patterns without team agreement.
