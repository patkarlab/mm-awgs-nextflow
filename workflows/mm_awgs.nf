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
}
