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
include { ICHORKARYO           } from '../../modules/local/ichorkaryo.nf'
include { GENEBAF              } from '../../modules/local/genebaf.nf'
include { VEP_ANNOTATE_CLAIR3  } from '../../modules/local/vep_annotate_clair3.nf'
include { FILTER_V6_REPORT    } from '../../modules/local/filter_v6_report.nf'
include { BAF_LOH_SCREEN      } from '../../modules/local/baf_loh_screen.nf'
include { BAF_CN_PLOTS        } from '../../modules/local/baf_cn_plots.nf'

workflow HG38_TRACK {

    take:
    minknow_bams   // [meta, minknow_bam]
    ch_panel_bed   // [meta, panel_bed]  hg38 panel this sample was sequenced with

    main:
    // REALIGN_HG38 emits [meta, bam, bai] directly — indexing is done in the
    // same work dir as the realignment to keep BAM and BAI colocated. This
    // matters for tools like Clair3 whose wrapper resolves the BAM symlink
    // and looks for the BAI next to the resolved path.
    REALIGN_HG38(minknow_bams)
    hg38_bam_bai = REALIGN_HG38.out.bam_bai   // [meta, bam, bai]

    if (!params.skip_clairs_to) {
        // Joined on meta so calling is restricted to the panel this sample was
        // sequenced with. Calling a v4 sample over the v7 BED would report
        // regions that were never enriched, at off-target depth.
        CLAIRS_TO(
            hg38_bam_bai
                .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .join(ch_panel_bed.map { meta, bed -> tuple(meta.id, bed) })
                .map { id, meta, bam, bai, bed -> tuple(meta, bam, bai, bed) }
        )
    }
    if (!params.skip_ichorcna) {
        ICHORCNA(hg38_bam_bai)

        // Depth-based karyotype from the ichorCNA segments. Joined on meta.id,
        // not on the whole meta map: a map join goes empty the moment meta
        // gains a field between a cached task and a live one, and the run then
        // reports success having produced nothing.
        if (!params.skip_ichorkaryo) {
            ICHORKARYO(
                // cna_seg and wig join with remainder: true for the same
                // reason params_file does. ICHORCNA runs with errorStrategy
                // 'ignore', so a sample can already disappear from this
                // channel by failing upstream; a plain join on an optional
                // output would add a second, quieter way for one to vanish.
                // Missing here means null qc fields, not a missing sample.
                ICHORCNA.out.seg.map { meta, seg -> tuple(meta.id, meta, seg) }
                    .join(ICHORCNA.out.params_file.map { meta, p -> tuple(meta.id, p) },
                          remainder: true)
                    .map { id, meta, seg, p -> tuple(id, meta, seg, p ?: []) }
                    .join(ICHORCNA.out.cna_seg.map { meta, c -> tuple(meta.id, c) },
                          remainder: true)
                    .map { id, meta, seg, p, c -> tuple(id, meta, seg, p, c ?: []) }
                    .join(ICHORCNA.out.wig.map { meta, w -> tuple(meta.id, w) },
                          remainder: true)
                    .map { id, meta, seg, p, c, w -> tuple(id, meta, seg, p, c, w ?: []) }
                    .join(ch_panel_bed.map { meta, bed -> tuple(meta.id, bed) })
                    .map { _id, meta, seg, p, c, w, bed ->
                           tuple(meta, seg, p, c, w, bed) }
            )
        }
    }

    // Allelic state per panel gene, from on-target reads. Independent of
    // ichorCNA: it reads the BAM directly and is not gated on skip_ichorcna.
    if (!params.skip_genebaf) {
        GENEBAF(
            hg38_bam_bai.map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .join(ch_panel_bed.map { meta, bed -> tuple(meta.id, bed) })
                .map { _id, meta, bam, bai, bed -> tuple(meta, bam, bai, bed) }
        )
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

        // Cohort-level: one BED for every sample, and it must be one where
        // every sample has the same capture status at every region. See
        // params.cohort_bed_hg38 in nextflow.config.
        BAF_LOH_SCREEN(
            ids_ch,
            clair3_ch,
            ichor_ch,
            file(params.cohort_bed_hg38, checkIfExists: true)
        )
        baf_screen_ch = BAF_LOH_SCREEN.out.screen

        if (!params.skip_baf_cn_plots) {
            BAF_CN_PLOTS(
                BAF_LOH_SCREEN.out.screen,
                BAF_LOH_SCREEN.out.sample_map,
                clair3_ch,
                ichor_ch,
                file(params.cohort_bed_hg38, checkIfExists: true)
            )
        }
    }

    emit:
    hg38_bam_bai             = hg38_bam_bai
    clairs_to_outdir         = params.skip_clairs_to    ? Channel.empty() : CLAIRS_TO.out.outdir
    ichorcna_outdir          = params.skip_ichorcna     ? Channel.empty() : ICHORCNA.out.outdir
    ichorkaryo_json          = (params.skip_ichorcna || params.skip_ichorkaryo)
                               ? Channel.empty() : ICHORKARYO.out.karyotype
    genebaf_genes            = params.skip_genebaf      ? Channel.empty() : GENEBAF.out.genes
    genebaf_summary          = params.skip_genebaf      ? Channel.empty() : GENEBAF.out.summary
    genebaf_armloh           = params.skip_genebaf      ? Channel.empty() : GENEBAF.out.armloh
    clair3_phased_outdir     = params.skip_clair3_phased ? Channel.empty() : CLAIR3_PHASED.out.outdir
    clair3_annotated_outdir  = (params.skip_clair3_phased || params.skip_vep_annotate) ? Channel.empty() : VEP_ANNOTATE_CLAIR3.out.outdir
    v6_report                = (params.skip_clair3_phased || params.skip_vep_annotate || params.skip_v6_filter) ? Channel.empty() : FILTER_V6_REPORT.out.clinical
    baf_loh_screen           = baf_screen_ch
}
