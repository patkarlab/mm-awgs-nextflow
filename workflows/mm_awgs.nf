/*
 * MM_AWGS top workflow.
 *
 * Dual-reference pipeline:
 *   Input MinKNOW BAM ─┬─► T2T track  (SV)
 *                      └─► hg38 track (SNV/CNV/germline)
 *
 * Each track owns its own realignment, so there is no inter-track dependency
 * and the two tracks run concurrently on Nextflow's scheduler.
 */

include { PREPARE_INPUT } from '../subworkflows/local/prepare_input.nf'
include { MERGE_MINKNOW } from '../modules/local/merge_minknow.nf'
include { T2T_TRACK     } from '../subworkflows/local/t2t_track.nf'
include { HG38_TRACK    } from '../subworkflows/local/hg38_track.nf'
include { REPORT_BUNDLE } from '../modules/local/report_bundle.nf'
include { DASHBOARD     } from '../modules/local/dashboard.nf'

workflow MM_AWGS {

    // 1. Sample sheet → channel of [meta, minknow_input]
    //    minknow_input is a FOLDER of per-chunk BAMs from the P2i (one BAM
    //    per ~1.5 h), or a single BAM (validation sheets).
    PREPARE_INPUT()

    // 2. Concatenate each sample's chunk BAMs into a single BAM. Done once per
    //    sample so both tracks reuse the same merged input.
    MERGE_MINKNOW(PREPARE_INPUT.out.minknow_bams)
    minknow_bams = MERGE_MINKNOW.out.merged_bam   // [meta, merged_bam]

    // 3. T2T track (SV calling + MM annotation)
    if (!params.skip_t2t_track) {
        T2T_TRACK(minknow_bams)
    }

    // 4. hg38 track (SNV/indel via ClairS-TO + large CNV via ichorCNA +
    //               phased germline via Clair3 + VEP annotation)
    if (!params.skip_hg38_track) {
        HG38_TRACK(minknow_bams)
    }

    // 5. Reporting.
    //
    // The bundle scans the published output tree, so it must not start until
    // publishDir has run for the processes that produce what it collects. The
    // outputs of both tracks are collected into a single signal channel used
    // only for ordering; REPORT_BUNDLE ignores its contents.
    if (!params.skip_report_bundle) {
        // Every emit the bundle collects from is mixed in, not just one. The
        // bundle copies SV annotations, translocations, SNV tables, ichorCNA
        // output, on-target QC and the BAF/LOH screen; waiting on only one of
        // those would let it start while the others were still publishing.
        ready = Channel.empty()
        if (!params.skip_t2t_track) {
            ready = ready
                .mix(T2T_TRACK.out.mm_annotated_tsv.map { it -> 'ok' })
                .mix(T2T_TRACK.out.translocations.map   { it -> 'ok' })
                .mix(T2T_TRACK.out.qc_coverage.map      { it -> 'ok' })
        }
        if (!params.skip_hg38_track) {
            ready = ready
                .mix(HG38_TRACK.out.ichorcna_outdir.map { it -> 'ok' })
                .mix(HG38_TRACK.out.v6_report.map       { it -> 'ok' })
                .mix(HG38_TRACK.out.baf_loh_screen.map  { it -> 'ok' })
        }

        bundle_name = params.report_bundle_name
            ?: "report_" + file(params.outdir).getName()

        REPORT_BUNDLE(
            ready.collect().ifEmpty(['none']),
            file("${projectDir}/bin/build_report_bundle.sh"),
            params.outdir,
            bundle_name
        )

        if (!params.skip_dashboard) {
            DASHBOARD(
                REPORT_BUNDLE.out.bundle,
                file("${projectDir}/bin/dashboard_builder")
            )
        }
    }
}
