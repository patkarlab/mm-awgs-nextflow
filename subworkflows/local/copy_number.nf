// RETIRED 2026-09-16 (retire_mmkaryo_v1). No longer included by any
// workflow. Copy-number figures come from ICHORCNA_PLOT on the hg38 track;
// see modules/local/ichorcna_plot.nf. Kept for one release so the history
// of the copy-number tab stays readable, then delete.
//
// COPY_NUMBER: allele-specific copy number from the off-target reads of the
// adaptive sampling run.
//
// Runs on the T2T track, because the panel BED, the centromere table and the
// cytoband file are all in CHM13v2.0 coordinates, and because T2T resolves the
// IGH locus that hg38 sequesters into alt contigs.
//
// This is a distinct capability from ichorCNA on the hg38 track. ichorCNA is a
// tumour-fraction estimator that emits copy number; this produces segmented
// copy number with allele-specific states, per-segment clonal cell fractions
// comparable to FISH percentages, and an explicit detection limit.
//
// Note what is deliberately NOT surfaced downstream: the baseline
// disomic/trisomic verdict and any free-purity fit. Both are documented as
// unreliable until BAF gets its own segmentation, and in a report they would
// read as authoritative.
//
// The panel BED arrives per sample rather than as a pipeline-wide value. The
// cohort spans four panel versions, and the BED is what mmkaryo masks: masking
// a v4 sample with the v7 panel discards genuine off-target signal and leaves
// v4-only enriched regions unmasked, which shifts the profile in both
// directions at once.
//
// Purity reaches the processes through meta.purity, set by PREPARE_INPUT from
// the sample sheet. A sample with no purity runs a free purity/ploidy grid fit,
// which is degenerate between ploidy basins, and PREPARE_INPUT warns by name
// when that happens.
//

include { MMKARYO } from '../../modules/local/mmkaryo'
include { MMBAF   } from '../../modules/local/mmbaf'
include { MMPLOT  } from '../../modules/local/mmplot'

workflow COPY_NUMBER {

    take:
    ch_t2t_bam        // channel: [ val(meta), path(bam), path(bai) ]
    ch_panel_bed      // channel: [ val(meta), path(panel_bed) ]  per-sample T2T panel
    t2t_fai           // path
    t2t_fasta         // path
    cytobands         // path (or NO_FILE)
    genes             // path (or NO_FILE)
    g1000_vcf         // path

    main:
    ch_versions = Channel.empty()

    // Join on meta so a sample can never be masked with another sample's panel.
    // join() is a strict inner join here on purpose: a sample present in one
    // channel and not the other is a sample sheet error, and dropping it
    // silently would be worse than the pipeline finishing short.
    // Keyed on meta.id, not on the whole meta map. join() defaults to matching
    // element 0, and when that is the map itself the join succeeds only while
    // every field agrees on both sides. A mismatch produces an empty channel,
    // which Nextflow treats as success: the processes simply never run.
    ch_mmkaryo_in = ch_t2t_bam
        .map { meta, bam, bai -> tuple(meta.id, meta, bam, bai) }
        .join(ch_panel_bed.map { meta, bed -> tuple(meta.id, bed) })
        .map { id, meta, bam, bai, bed -> tuple(meta, bam, bai, bed) }

    MMKARYO(ch_mmkaryo_in, t2t_fai, t2t_fasta, cytobands)
    ch_versions = ch_versions.mix(MMKARYO.out.versions.first())

    // BAF needs the BAM again alongside mmkaryo's segments and fit.
    ch_baf_in = ch_t2t_bam
        .join(MMKARYO.out.segments)
        .join(MMKARYO.out.karyotype)

    MMBAF(ch_baf_in, t2t_fai, g1000_vcf)
    ch_versions = ch_versions.mix(MMBAF.out.versions.first())

    // The first three joins pair outputs of processes that received the same
    // meta object, so matching on the map is safe there. The BED comes from
    // PREPARE_INPUT and is keyed on meta.id for the reason above.
    ch_plot_in = MMKARYO.out.bins
        .join(MMBAF.out.segments)
        .join(MMBAF.out.windows)
        .join(MMKARYO.out.karyotype)
        .map { meta, bins, baf, win, karyo -> tuple(meta.id, meta, bins, baf, win, karyo) }
        .join(ch_panel_bed.map { meta, bed -> tuple(meta.id, bed) })
        .map { id, meta, bins, baf, win, karyo, bed ->
            tuple(meta, bins, baf, win, karyo, bed)
        }

    MMPLOT(ch_plot_in, cytobands, genes)
    ch_versions = ch_versions.mix(MMPLOT.out.versions.first())

    emit:
    bins        = MMKARYO.out.bins
    segments    = MMKARYO.out.segments
    arms        = MMKARYO.out.arms
    karyotype   = MMKARYO.out.karyotype
    cytobands   = MMKARYO.out.cytobands
    baf         = MMBAF.out.segments
    baf_windows = MMBAF.out.windows
    baf_summary = MMBAF.out.summary
    plot_genome = MMPLOT.out.genome
    plot_grid   = MMPLOT.out.grid
    versions    = ch_versions
}
