/*
 * PREPARE_INPUT
 *
 * Reads the sample sheet CSV and emits a channel of [meta, minknow_input].
 *
 * Expected sample sheet columns (header row required):
 *   sample_id      Canonical sequencing ID (e.g. 11F20262905). No PHI.
 *   minknow_bam    Path to the sample's P2i output. This may be either:
 *                    - a FOLDER of per-chunk MinKNOW BAMs (real P2i output;
 *                      one BAM per ~1.5 h, NC_-named, T2T-aligned), OR
 *                    - a single BAM file (validation sheets).
 *                  MERGE_MINKNOW concatenates a folder into one BAM; a single
 *                  file passes through. Point a folder column at the directory
 *                  that directly contains the pass BAMs (e.g. .../bam_pass).
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
            // checkIfExists: true so a bad path (file OR directory) fails fast
            // and legibly at launch, rather than deep inside MERGE_MINKNOW.
            tuple(meta, file(row.minknow_bam.trim(), checkIfExists: true))
        }

    emit:
    minknow_bams = minknow_bams
}
