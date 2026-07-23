/*
 * HG38_TRACK
 *
 * Realign MinKNOW BAM to hg38, then in parallel:
 *   - ClairS-TO  : somatic SNV/indel (docker)
 *   - ichorCNA   : large CNV (native conda env, per-sample)
 *   - Clair3     : phased germline VCF (docker)
 *
 * Clair3 output feeds VEP for annotation.
 *
 * Input:  channel [meta, minknow_bam]
 * Output:
 *   hg38_bam_bai           [meta, bam, bai]
 *   clairs_to_outdir       [meta, outdir]
 *   ichorcna_outdir        [meta, outdir]
 *   clair3_phased_outdir   [meta, outdir]
 *   clair3_annotated_outdir [meta, outdir]
 */

include { REALIGN_HG38         } from '../../modules/local/realign_hg38.nf'
include { CLAIRS_TO            } from '../../modules/local/clairs_to.nf'
include { ICHORCNA             } from '../../modules/local/ichorcna.nf'
include { CLAIR3_PHASED        } from '../../modules/local/clair3_phased.nf'
include { VEP_ANNOTATE_CLAIR3  } from '../../modules/local/vep_annotate_clair3.nf'
include { FILTER_V6_REPORT    } from '../../modules/local/filter_v6_report.nf'
include { BAF_LOH_SCREEN      } from '../../modules/local/baf_loh_screen.nf'
include { BAF_CN_PLOTS        } from '../../modules/local/baf_cn_plots.nf'

workflow HG38_TRACK {

    take:
    minknow_bams   // [meta, minknow_bam]

    main:
    // REALIGN_HG38 emits [meta, bam, bai] directly — indexing is done in the
    // same work dir as the realignment to keep BAM and BAI colocated. This
    // matters for tools like Clair3 whose wrapper resolves the BAM symlink
    // and looks for the BAI next to the resolved path.
    REALIGN_HG38(minknow_bams)
    hg38_bam_bai = REALIGN_HG38.out.bam_bai   // [meta, bam, bai]

    if (!params.skip_clairs_to) {
        CLAIRS_TO(hg38_bam_bai)
    }
    if (!params.skip_ichorcna) {
        ICHORCNA(hg38_bam_bai)
    }
    if (!params.skip_clair3_phased) {
        CLAIR3_PHASED(hg38_bam_bai)
        if (!params.skip_vep_annotate) {
            VEP_ANNOTATE_CLAIR3(CLAIR3_PHASED.out.merge_output)
            if (!params.skip_v6_filter) {
                FILTER_V6_REPORT(VEP_ANNOTATE_CLAIR3.out.candidates_tsv)
            }
        }
    }

    // Cohort BAF / LOH screen. Collected rather than per-sample: the screen
    // compares each region against the cohort median for that same region, and
    // that baseline is withheld below three samples.
    //
    // Sample identifiers, phased VCF directories and ichorCNA directories are
    // collected separately and matched inside the process by inspecting the
    // staged files, since Nextflow does not guarantee that separately collected
    // channels stage in the same order.
    baf_screen_ch = Channel.empty()
    if (!params.skip_clair3_phased && !params.skip_baf_loh) {
        ids_ch    = CLAIR3_PHASED.out.outdir.map { meta, _dir -> meta.id }.collect()
        clair3_ch = CLAIR3_PHASED.out.outdir.map { _meta, dir -> dir }.collect()
        ichor_ch  = params.skip_ichorcna ?
            Channel.value([]) :
            ICHORCNA.out.outdir.map { _meta, dir -> dir }.collect().ifEmpty([])

        BAF_LOH_SCREEN(
            ids_ch,
            clair3_ch,
            ichor_ch,
            file(params.panel_bed_hg38)
        )
        baf_screen_ch = BAF_LOH_SCREEN.out.screen

        if (!params.skip_baf_cn_plots) {
            BAF_CN_PLOTS(
                BAF_LOH_SCREEN.out.screen,
                BAF_LOH_SCREEN.out.sample_map,
                clair3_ch,
                ichor_ch,
                file(params.panel_bed_hg38)
            )
        }
    }

    emit:
    hg38_bam_bai             = hg38_bam_bai
    clairs_to_outdir         = params.skip_clairs_to    ? Channel.empty() : CLAIRS_TO.out.outdir
    ichorcna_outdir          = params.skip_ichorcna     ? Channel.empty() : ICHORCNA.out.outdir
    clair3_phased_outdir     = params.skip_clair3_phased ? Channel.empty() : CLAIR3_PHASED.out.outdir
    clair3_annotated_outdir  = (params.skip_clair3_phased || params.skip_vep_annotate) ? Channel.empty() : VEP_ANNOTATE_CLAIR3.out.outdir
    v6_report                = (params.skip_clair3_phased || params.skip_vep_annotate || params.skip_v6_filter) ? Channel.empty() : FILTER_V6_REPORT.out.clinical
    baf_loh_screen           = baf_screen_ch
}
