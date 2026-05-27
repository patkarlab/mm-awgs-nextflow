# Usage

## Sample sheet

CSV with header row and these columns:

| Column        | Required | Description |
|---------------|----------|-------------|
| `sample_id`   | yes      | Canonical sequencing ID. No PHI. |
| `minknow_bam` | yes      | Path to the raw MinKNOW BAM (NC_-named, T2T-aligned at P2i). |
| `timepoint`   | no       | Free text. `18h` is recognized as an early snapshot. Other values → `full`. |
| `notes`       | no       | Free text. Ignored by the pipeline. |

Template: [`assets/sample_sheet_template.csv`](../assets/sample_sheet_template.csv).

## Required runtime input

Just `--sample_sheet`. All reference paths default to gandalf locations
(see `nextflow.config` parameter defaults). Override any when running off-server.

## Profiles

Profiles compose. The recommended combo on the project server is:

```bash
-profile conda,docker,gandalf
```

| Profile     | Effect |
|-------------|--------|
| `conda`     | Enable conda for any process with a `conda` directive. |
| `docker`    | Enable docker for any process with a `container` directive. |
| `gandalf`   | Local executor, 64 CPU / 128 GB, work dir on `/goast`. |
| `test`      | Tiny synthetic inputs. |
| `stub`      | Stub-run; each module's `stub:` block runs instead of real tool. |
| `singularity` | Use singularity instead of docker (untested in v0.1). |

## Track and step toggles

The pipeline has many `--skip_*` flags so you can iterate on one component:

| Flag                       | Effect |
|----------------------------|--------|
| `--skip_t2t_track`         | Skip the entire T2T (SV) track |
| `--skip_hg38_track`        | Skip the entire hg38 (SNV/CNV) track |
| `--skip_sv_calling`        | Within T2T: skip callers + SURVIVOR + annotation |
| `--skip_mm_annotation`     | Within T2T: skip only MM annotation, keep SV calling |
| `--skip_clairs_to`         | Within hg38: skip ClairS-TO |
| `--skip_ichorcna`          | Within hg38: skip ichorCNA |
| `--skip_clair3_phased`     | Within hg38: skip Clair3 (and VEP cascades off) |
| `--skip_vep_annotate`      | Within hg38: skip VEP only; Clair3 still runs |

## Tuning knobs

The default values match the v3-era production tuning. Override any:

```bash
# E.g. raise minimum caller support for less noisy SV calls
nextflow run main.nf -profile conda,docker,gandalf \
    --sample_sheet samples.csv \
    --sniffles_min_support 4 \
    --cutesv_min_support 4 \
    --severus_min_support 3
```

All knobs are listed in `nextflow.config` under `params { ... }`. The main
groups:

- Reference paths: `t2t_fasta`, `t2t_mmi`, `hg38_fasta`, `hg38_mmi`
- Panel BEDs: `panel_bed_t2t`, `panel_bed_hg38`
- Sniffles: `sniffles_min_support`, `sniffles_min_mapq`, `sniffles_min_aln_len`
- CuteSV: `cutesv_min_support`, `cutesv_min_mapq`, `cutesv_min_size`, `cutesv_min_read_len`
- Severus: `severus_min_support`, `severus_min_mapq`, `severus_vntr_bed`, `severus_pon`
- SURVIVOR: `survivor_max_dist`, `survivor_min_callers`, `survivor_take_type`, `survivor_take_strand`, `survivor_estimate_dist`, `survivor_min_size`
- ClairS-TO: `clairs_to_platform`, `clairs_to_threads`, `clairs_to_disable_verdict`
- ichorCNA: `ichorcna_bin_size`, `ichorcna_ploidy`, `ichorcna_normal`, `ichorcna_max_cn`, `ichorcna_chrtrain`
- Clair3: `clair3_model_name`, `clair3_threads`, `clair3_platform`
- VEP: `vep_cache_version`, `vep_gnomad_af_max`, `vep_exclude_consequences`

## Resume

`-resume` works out of the box. Cached tasks are reused unless their inputs
or parameters change.

## Output directory

Each track publishes under `${params.outdir}/{t2t,hg38}/...`. See
[`output.md`](output.md) for the full layout.
