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
include { SAMTOOLS_INDEX_HG38  } from '../../modules/local/samtools_index.nf'
include { CLAIRS_TO            } from '../../modules/local/clairs_to.nf'
include { ICHORCNA             } from '../../modules/local/ichorcna.nf'
include { CLAIR3_PHASED        } from '../../modules/local/clair3_phased.nf'
include { VEP_ANNOTATE_CLAIR3  } from '../../modules/local/vep_annotate_clair3.nf'

workflow HG38_TRACK {

    take:
    minknow_bams   // [meta, minknow_bam]

    main:
    REALIGN_HG38(minknow_bams)
    SAMTOOLS_INDEX_HG38(REALIGN_HG38.out.bam)
    hg38_bam_bai = SAMTOOLS_INDEX_HG38.out.bam_bai   // [meta, bam, bai]

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
        }
    }

    emit:
    hg38_bam_bai             = hg38_bam_bai
    clairs_to_outdir         = params.skip_clairs_to    ? Channel.empty() : CLAIRS_TO.out.outdir
    ichorcna_outdir          = params.skip_ichorcna     ? Channel.empty() : ICHORCNA.out.outdir
    clair3_phased_outdir     = params.skip_clair3_phased ? Channel.empty() : CLAIR3_PHASED.out.outdir
    clair3_annotated_outdir  = (params.skip_clair3_phased || params.skip_vep_annotate) ? Channel.empty() : VEP_ANNOTATE_CLAIR3.out.outdir
}
