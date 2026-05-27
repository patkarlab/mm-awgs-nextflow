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
include { T2T_TRACK     } from '../subworkflows/local/t2t_track.nf'
include { HG38_TRACK    } from '../subworkflows/local/hg38_track.nf'

workflow MM_AWGS {

    // 1. Sample sheet → channel of [meta, minknow_bam]
    PREPARE_INPUT()
    minknow_bams = PREPARE_INPUT.out.minknow_bams

    // 2. T2T track (SV calling + MM annotation)
    if (!params.skip_t2t_track) {
        T2T_TRACK(minknow_bams)
    }

    // 3. hg38 track (SNV/indel via ClairS-TO + large CNV via ichorCNA +
    //               phased germline via Clair3 + VEP annotation)
    if (!params.skip_hg38_track) {
        HG38_TRACK(minknow_bams)
    }
}
