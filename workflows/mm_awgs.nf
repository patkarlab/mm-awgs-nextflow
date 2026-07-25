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
include { IGV_SNAPSHOTS       } from '../modules/local/igv_snapshots.nf'
include { EMBED_REPORT_ASSETS } from '../modules/local/embed_report_assets.nf'
include { REPORT_ZIP          } from '../modules/local/report_zip.nf'

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

    // 4b. IGV snapshots.
    //
    // Two evidence classes per sample: paired breakpoint pages for each
    // rearrangement against T2T, and one clinical SNV page against hg38. The
    // SNV page is published under the filename the dashboard builder
    // resolves, which also gives the variant cards their IGV links.
    //
    // Joined on meta so each sample's tables stay with its own alignments.
    // v6_report uses remainder because a sample with no on-panel clinical
    // SNVs never emits one; without it that sample would be dropped here and
    // lose its translocation pages too.
    if (!params.skip_igv && !params.skip_t2t_track && !params.skip_hg38_track) {
        igv_input = T2T_TRACK.out.mm_annotated_tsv
            .join(T2T_TRACK.out.t2t_bam_bai)
            .join(HG38_TRACK.out.hg38_bam_bai)
            .join(HG38_TRACK.out.v6_report, remainder: true)
            .filter { it[0] != null && it[1] != null }
            .map { meta, mm, tbam, tbai, hbam, hbai, clin ->
                tuple(meta, mm, clin ?: [], tbam, tbai, hbam, hbai)
            }

        IGV_SNAPSHOTS(
            igv_input,
            file(params.t2t_fasta),
            file(params.t2t_fai),
            file(params.hg38_fasta),
            file(params.hg38_fai)
        )
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

        // IGV pages are collected by the bundle, so the bundle must not
        // start before they are published.
        if (!params.skip_igv && !params.skip_t2t_track && !params.skip_hg38_track) {
            ready = ready.mix(IGV_SNAPSHOTS.out.igv.map { it -> 'ok' })
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

            // Inline every local dependency, then package. Without the embed
            // step the reports reference their stylesheets and figures by
            // relative path and break as soon as they are moved.
            if (!params.skip_report_package) {
                EMBED_REPORT_ASSETS(DASHBOARD.out.bundle)
                REPORT_ZIP(EMBED_REPORT_ASSETS.out.bundle)
            }
        }
    }
}
