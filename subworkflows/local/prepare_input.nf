/*
 * PREPARE_INPUT
 *
 * Reads the sample sheet CSV and emits a channel of [meta, minknow_bam].
 *
 * Expected sample sheet columns (header row required):
 *   sample_id      Canonical sequencing ID (e.g. 11F20262905). No PHI.
 *   minknow_bam    Path to the raw MinKNOW BAM (NC_-named, T2T-aligned at P2i).
 *   timepoint      Optional. Free text; '18h' is recognized as a paired snapshot
 *                  marker by downstream cohort QC. Other values treated as 'full'.
 *   notes          Optional. Free text. Ignored by the pipeline.
 */

workflow PREPARE_INPUT {

    main:
    if (!params.sample_sheet) {
        error "sample_sheet not set; provide --sample_sheet path/to/samples.csv"
    }

    minknow_bams = Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            if (!row.sample_id) {
                error "sample sheet row missing sample_id: ${row}"
            }
            if (!row.minknow_bam) {
                error "sample sheet row missing minknow_bam for ${row.sample_id}"
            }
            def meta = [
                id        : row.sample_id.trim(),
                timepoint : (row.timepoint ?: '').trim() ?: null,
                notes     : (row.notes ?: '').trim() ?: null,
            ]
            tuple(meta, file(row.minknow_bam.trim(), checkIfExists: false))
        }

    emit:
    minknow_bams = minknow_bams
}
