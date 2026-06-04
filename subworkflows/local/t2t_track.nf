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

workflow T2T_TRACK {

    take:
    minknow_bams   // [meta, minknow_bam]

    main:
    // REALIGN_T2T emits [meta, bam, bai] directly — indexing is done in the
    // same work dir as the realignment to keep BAM and BAI colocated for
    // downstream tools that resolve symlinks (e.g. Clair3's wrapper).
    REALIGN_T2T(minknow_bams)
    t2t_bam_bai = REALIGN_T2T.out.bam_bai   // [meta, bam, bai]

    if (!params.skip_sv_calling) {
        SNIFFLES(t2t_bam_bai)
        CUTESV(t2t_bam_bai)
        SEVERUS(t2t_bam_bai)

        // Combine per-sample VCFs for merging
        per_sample_for_merge = SNIFFLES.out.vcf
            .join(CUTESV.out.vcf,      by: 0)
            .join(SEVERUS.out.vcf,     by: 0)
            .map { meta, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf ->
                tuple(meta, sn_vcf, cu_vcf, sv_vcf)
            }

        SURVIVOR_MERGE(per_sample_for_merge)

        if (!params.skip_mm_annotation) {
            ANNOTATE_MM_TRANSLOCATIONS(SURVIVOR_MERGE.out.merged_vcf)

            // Layer real per-caller read support onto the annotated TSV,
            // then unite near-identical translocation calls. Join on meta;
            // drop the .tbi from Sniffles/CuteSV (augmenter needs the .vcf);
            // Severus is already a plain .vcf.
            ch_augment_in = ANNOTATE_MM_TRANSLOCATIONS.out.tsv
                .join(SNIFFLES.out.vcf, by: 0)
                .join(CUTESV.out.vcf,   by: 0)
                .join(SEVERUS.out.vcf,  by: 0)
                .map { meta, annotated, sn_vcf, sn_tbi, cu_vcf, cu_tbi, sv_vcf ->
                    tuple(meta, annotated, sn_vcf, cu_vcf, sv_vcf)
                }

            AUGMENT_SV_SUPPORT(ch_augment_in)
            MERGE_TRANSLOCATIONS(AUGMENT_SV_SUPPORT.out.annotated)
        }
    }

    emit:
    t2t_bam_bai      = t2t_bam_bai
    sniffles_vcf     = params.skip_sv_calling     ? Channel.empty() : SNIFFLES.out.vcf
    cutesv_vcf       = params.skip_sv_calling     ? Channel.empty() : CUTESV.out.vcf
    severus_outdir   = params.skip_sv_calling     ? Channel.empty() : SEVERUS.out.outdir
    merged_vcf       = params.skip_sv_calling     ? Channel.empty() : SURVIVOR_MERGE.out.merged_vcf
    mm_annotated_tsv = (params.skip_sv_calling || params.skip_mm_annotation) ? Channel.empty() : AUGMENT_SV_SUPPORT.out.annotated
    translocations   = (params.skip_sv_calling || params.skip_mm_annotation) ? Channel.empty() : MERGE_TRANSLOCATIONS.out.translocations
}
