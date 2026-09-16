/*
 * T2T_TRACK
 *
 * Realign MinKNOW BAM to T2T-CHM13v2.0, then run three SV callers in
 * parallel, merge with SURVIVOR, and annotate against the MM partner-pair
 * dictionary.
 *
 * Input:  channel [meta, minknow_bam]
 * Output:
 *   t2t_bam_bai           [meta, bam, bai]
 *   sniffles_vcf          [meta, vcf, tbi]
 *   cutesv_vcf            [meta, vcf, tbi]
 *   severus_outdir        [meta, outdir]
 *   merged_vcf            [meta, merged.vcf.gz, tbi]
 *   mm_annotated_tsv      [meta, tsv]
 */

include { REALIGN_T2T                 } from '../../modules/local/realign_t2t.nf'
include { SNIFFLES                    } from '../../modules/local/sniffles.nf'
include { CUTESV                      } from '../../modules/local/cutesv.nf'
include { SEVERUS                     } from '../../modules/local/severus.nf'
include { SURVIVOR_MERGE              } from '../../modules/local/survivor_merge.nf'
include { ANNOTATE_MM_TRANSLOCATIONS  } from '../../modules/local/annotate_mm_translocations.nf'
include { AUGMENT_SV_SUPPORT   } from '../../modules/local/augment_sv_support'
include { MERGE_TRANSLOCATIONS } from '../../modules/local/merge_translocations'
include { QC_ONTARGET          } from '../../modules/local/qc_ontarget'
include { SAVANA               } from '../../modules/local/savana'
include { NANOMONSV            } from '../../modules/local/nanomonsv'   // nanomonsv_confirmation_v1

workflow T2T_TRACK {

    take:
    minknow_bams   // [meta, minknow_bam]
    ch_panel_bed   // [meta, panel_bed]  T2T panel this sample was sequenced with

    main:
    // REALIGN_T2T emits [meta, bam, bai] directly — indexing is done in the
    // same work dir as the realignment to keep BAM and BAI colocated for
    // downstream tools that resolve symlinks (e.g. Clair3's wrapper).
    REALIGN_T2T(minknow_bams)
    t2t_bam_bai = REALIGN_T2T.out.bam_bai   // [meta, bam, bai]

    // On-target QC: per-region coverage, read-length + basecaller-Q stats,
    // and histogram PNGs. Independent of SV calling; runs off the realigned
    // BAM and the v6 panel BED. Gated so it can be skipped.
    if (!params.skip_qc) {
        // Joined on meta: on-target coverage is only meaningful against the
        // panel this sample was actually enriched with. Measuring a v4 sample
        // against the v7 BED is what made every non-v7 coverage figure in the
        // cohort uninterpretable.
        QC_ONTARGET(
            t2t_bam_bai
                .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
                .join(ch_panel_bed.map { meta, bed -> tuple(meta.id, bed) })
                .map { id, meta, bam, bai, bed -> tuple(meta, bam, bai, bed) }
        )
    }

    if (!params.skip_sv_calling) {
        SNIFFLES(t2t_bam_bai)
        CUTESV(t2t_bam_bai)
        SEVERUS(t2t_bam_bai)
        SAVANA(t2t_bam_bai)

        // Confirmation layer, not a fifth vote: nanomonsv runs beside the
        // ensemble and its junctions are matched onto the annotated table
        // downstream. It never enters the SURVIVOR merge, so SUPP_VEC and
        // the callers column keep their meaning.
        nanomonsv_ch = Channel.empty()
        if (!params.skip_nanomonsv) {
            NANOMONSV(t2t_bam_bai)
            nanomonsv_ch = NANOMONSV.out.result
        }

        // Combine per-sample VCFs for merging
        per_sample_for_merge = SNIFFLES.out.vcf
            .join(CUTESV.out.vcf,      by: 0)
            .join(SEVERUS.out.vcf,     by: 0)
            .join(SAVANA.out.vcf,      by: 0)
            .map { meta, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf, sa_vcf ->
                tuple(meta, sn_vcf, cu_vcf, sv_vcf, sa_vcf)
            }

        SURVIVOR_MERGE(per_sample_for_merge)

        if (!params.skip_mm_annotation) {
            ANNOTATE_MM_TRANSLOCATIONS(SURVIVOR_MERGE.out.merged_vcf)

            // Layer real per-caller read support onto the annotated TSV,
            // then unite near-identical translocation calls. Join on meta;
            // drop the .tbi from Sniffles/CuteSV (augmenter needs the .vcf);
            // Severus is already a plain .vcf.
            // The nanomonsv result is joined with remainder: true on meta.id,
            // so a sample nanomonsv fails on (or a run with it skipped)
            // still gets its annotated table, with the confirmation columns
            // blank. A plain join would drop the sample silently.
            ch_augment_in = ANNOTATE_MM_TRANSLOCATIONS.out.tsv
                .join(SNIFFLES.out.vcf, by: 0)
                .join(CUTESV.out.vcf,   by: 0)
                .join(SEVERUS.out.vcf,  by: 0)
                .join(SAVANA.out.vcf,   by: 0)
                .map { meta, annotated, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf, sa_vcf ->
                    tuple(meta.id, meta, annotated, sn_vcf, cu_vcf, sv_vcf, sa_vcf)
                }
                .join(nanomonsv_ch.map { meta, res -> tuple(meta.id, res) }, remainder: true)
                .filter { it[1] != null }
                .map { _id, meta, annotated, sn_vcf, cu_vcf, sv_vcf, sa_vcf, nano ->
                    tuple(meta, annotated, sn_vcf, cu_vcf, sv_vcf, sa_vcf, nano ?: [])
                }

            AUGMENT_SV_SUPPORT(ch_augment_in)
            MERGE_TRANSLOCATIONS(AUGMENT_SV_SUPPORT.out.annotated)
        }
    }

    emit:
    t2t_bam_bai      = t2t_bam_bai
    qc_coverage      = params.skip_qc ? Channel.empty() : QC_ONTARGET.out.coverage
    qc_summary       = params.skip_qc ? Channel.empty() : QC_ONTARGET.out.summary
    sniffles_vcf     = params.skip_sv_calling     ? Channel.empty() : SNIFFLES.out.vcf
    cutesv_vcf       = params.skip_sv_calling     ? Channel.empty() : CUTESV.out.vcf
    severus_outdir   = params.skip_sv_calling     ? Channel.empty() : SEVERUS.out.outdir
    // The unfiltered classified VCF. savana_somatic_vcf is emitted alongside
    // it for comparison and is deliberately not consumed by the ensemble.
    savana_vcf       = params.skip_sv_calling     ? Channel.empty() : SAVANA.out.vcf
    savana_somatic_vcf = params.skip_sv_calling   ? Channel.empty() : SAVANA.out.somatic_vcf
    savana_outdir    = params.skip_sv_calling     ? Channel.empty() : SAVANA.out.outdir
    merged_vcf       = params.skip_sv_calling     ? Channel.empty() : SURVIVOR_MERGE.out.merged_vcf
    nanomonsv_result = (params.skip_sv_calling || params.skip_nanomonsv) ? Channel.empty() : NANOMONSV.out.result
    mm_annotated_tsv = (params.skip_sv_calling || params.skip_mm_annotation) ? Channel.empty() : AUGMENT_SV_SUPPORT.out.annotated
    translocations   = (params.skip_sv_calling || params.skip_mm_annotation) ? Channel.empty() : MERGE_TRANSLOCATIONS.out.translocations
}
