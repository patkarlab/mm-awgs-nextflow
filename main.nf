#!/usr/bin/env nextflow
/*
 * mm-awgs-nextflow
 * ----------------
 * Dual-reference adaptive WGS pipeline for multiple myeloma.
 *
 * T2T track:  realign → Sniffles2 + CuteSV + Severus → SURVIVOR → MM annotation
 * hg38 track: realign → ClairS-TO (SNV/indel) + ichorCNA (large CNV) +
 *                       Clair3 phased germline → VEP annotation
 *
 * Both tracks consume the same MinKNOW BAM and run in parallel.
 */

nextflow.enable.dsl = 2

// -----------------------------------------------------------------------------
// Help and version
// -----------------------------------------------------------------------------
def printHelp() {
    log.info """
    mm-awgs-nextflow v${workflow.manifest.version}

    USAGE
        nextflow run main.nf \\
            -profile conda,docker,gandalf \\
            --sample_sheet samples.csv \\
            --outdir       results/

    REQUIRED
        --sample_sheet      CSV: sample_id,minknow_bam,timepoint,notes
                            (see assets/sample_sheet_template.csv)

    REFERENCES (defaults point to gandalf paths; override if running elsewhere)
        --t2t_fasta         T2T-CHM13v2.0 chr-named FASTA
        --t2t_mmi           Pre-built minimap2 index (optional but recommended)
        --hg38_fasta        hg38 FASTA
        --severus_vntr_bed  T2T VNTR BED
        --severus_pon       T2T population panel-of-normals

    PANEL BEDS (defaults from assets/)
        --panel_bed_t2t     v6 T2T panel (chr-named) for SV annotation
        --panel_bed_hg38    v5 hg38 panel for ClairS-TO calling region
                            [NOTE: v6 hg38 BED is a v0.2 task]

    TRACK TOGGLES
        --skip_t2t_track          Skip the entire T2T (SV) track
        --skip_hg38_track         Skip the entire hg38 (SNV/CNV) track
        --skip_sv_calling         Within T2T: skip Sniffles+CuteSV+Severus+merge
        --skip_mm_annotation      Within T2T: skip MM-specific SV annotation
        --skip_clairs_to          Within hg38: skip ClairS-TO SNV/indel calling
        --skip_ichorcna           Within hg38: skip ichorCNA large CNV
        --skip_clair3_phased      Within hg38: skip Clair3 phased germline
        --skip_vep_annotate       Within hg38: skip VEP annotation of Clair3

    PROFILES
        -profile conda,docker,gandalf       recommended on the project server
        -profile conda                      conda envs only (some steps need docker)
        -profile docker                     docker only
        -profile gandalf                    gandalf-specific (local executor, work on /goast)
        -profile test                       small synthetic inputs
        -profile stub                       stub-only mode (CI smoke test)

    OUTPUT
        --outdir            Output directory [${params.outdir}]
    """
}

if (params.help) { printHelp(); exit 0 }
if (params.version) {
    log.info "mm-awgs-nextflow v${workflow.manifest.version}"
    exit 0
}

// -----------------------------------------------------------------------------
// Input validation
// -----------------------------------------------------------------------------
def required = [sample_sheet: params.sample_sheet]
def missing = required.findAll { _k, v -> v == null || v.toString().trim().isEmpty() }
if (missing) {
    log.error "Missing required parameters: ${missing.keySet().join(', ')}"
    log.error "Run with --help for details."
    exit 1
}

// -----------------------------------------------------------------------------
// Workflow
// -----------------------------------------------------------------------------
include { MM_AWGS } from './workflows/mm_awgs.nf'

workflow {
    MM_AWGS()
}

workflow.onComplete {
    log.info """
    ${'='*60}
    mm-awgs-nextflow run complete
    ${'='*60}
    status:        ${workflow.success ? 'SUCCESS' : 'FAILED'}
    runName:       ${workflow.runName}
    duration:      ${workflow.duration}
    outdir:        ${params.outdir}
    workDir:       ${workflow.workDir}
    """.stripIndent()
}

workflow.onError {
    log.error "Pipeline failed: ${workflow.errorMessage}"
}
