/*
 * PREPARE_INPUT
 *
 * Reads the sample sheet CSV and emits:
 *   minknow_bams  [meta, minknow_input]
 *   panel_beds    [meta, panel_bed_t2t, panel_bed_hg38]
 *
 * Expected sample sheet columns (header row required):
 *   sample_id      Canonical sequencing ID. No PHI.
 *   minknow_bam    Path to the sample's P2i output. Either:
 *                    - a FOLDER of per-chunk MinKNOW BAMs (real P2i output;
 *                      one BAM per ~1.5 h, NC_-named, T2T-aligned), OR
 *                    - a single BAM file (validation sheets).
 *                  MERGE_MINKNOW concatenates a folder into one BAM; a single
 *                  file passes through. Point a folder column at the directory
 *                  that directly contains the pass BAMs.
 *   purity         Optional. Tumour purity as a FRACTION in (0, 1], from
 *                  post-sort flow cytometry. See the note below.
 *   panel_version  Optional. One of v4, v5, v6, v7. Selects the adaptive
 *                  sampling panel BED this sample was actually sequenced with.
 *                  Falls back to params.panel_bed_t2t / params.panel_bed_hg38.
 *   timepoint      Optional. Free text; '18h' is recognized as a paired snapshot
 *                  marker by downstream cohort QC. Other values treated as 'full'.
 *   notes          Optional. Free text. Ignored by the pipeline.
 *
 * Why purity is validated rather than passed through
 * --------------------------------------------------
 * mmkaryo and mmbaf both take purity as `meta.purity ? "--purity ..." : ''`.
 * A null therefore produces a free purity/ploidy grid fit with no error and no
 * log line. That fit is structurally degenerate: (purity 0.30, ploidy 2.0) and
 * (purity 0.43, ploidy 4.0) score nearly identically, so a free fit can return
 * the wrong ploidy basin and every integer copy state follows from it. Until
 * this sheet carried a purity column, the whole cohort ran free.
 *
 * mmkaryo.py clamps a supplied purity into [0.01, 0.99], so a value entered as
 * a percentage (93.27) would be silently pinned at 0.99. Purity is therefore
 * required to be a fraction here, and anything outside (0, 1] is a hard error
 * at launch rather than a wrong number in a report.
 *
 * A sample with no purity still runs, with a named warning. Its copy number is
 * not directly comparable to the pinned samples and should be reported as such.
 */

// Panel version to BED basename, for both references. The v4-v6 panels are the
// MM-focused series; v7 is the broader plasma-cell-neoplasm panel.
def panelBedNames = [
    v4: [t2t: 'aWGS_MMfocused_v4_t2t_chr.bed', hg38: 'aWGS_MMfocused_v4_hg38.bed'],
    v5: [t2t: 'aWGS_MMfocused_v5_t2t_chr.bed', hg38: 'aWGS_MMfocused_v5_hg38.bed'],
    v6: [t2t: 'aWGS_MMfocused_v6_t2t_chr.bed', hg38: 'aWGS_MMfocused_v6_hg38.bed'],
    v7: [t2t: 'aWGS_PCN_v7_t2t_chr.bed',       hg38: 'aWGS_PCN_v7_hg38.bed'],
]

workflow PREPARE_INPUT {

    main:
    if (!params.sample_sheet) {
        error "sample_sheet not set; provide --sample_sheet path/to/samples.csv"
    }

    ch_rows = Channel
        .fromPath(params.sample_sheet, checkIfExists: true)
        .splitCsv(header: true, sep: ',')
        .map { row ->

            if (!row.sample_id) {
                error "sample sheet row missing sample_id: ${row}"
            }
            def sample_id = row.sample_id.trim()

            if (!row.minknow_bam) {
                error "sample sheet row missing minknow_bam for ${sample_id}"
            }

            // --- purity -------------------------------------------------
            // Empty is allowed and means "not measured". A present value must
            // be a fraction; a percentage is an error, not something to coerce.
            def purity = null
            def purity_text = (row.purity ?: '').trim()
            if (purity_text) {
                if (!purity_text.isNumber()) {
                    error "sample ${sample_id}: purity '${purity_text}' is not numeric. " +
                          "Supply a fraction in (0, 1], or leave blank if not measured."
                }
                def purity_value = purity_text.toDouble()
                if (purity_value <= 0.0 || purity_value > 1.0) {
                    error "sample ${sample_id}: purity ${purity_value} is outside (0, 1]. " +
                          "Purity must be a FRACTION, not a percentage: 93.27 percent is 0.9327. " +
                          "mmkaryo clamps to [0.01, 0.99], so an unconverted percentage " +
                          "would be pinned at 0.99 without any error."
                }
                purity = purity_value
            }
            else {
                log.warn "sample ${sample_id}: no purity in the sample sheet. " +
                         "mmkaryo and mmbaf will run a FREE purity/ploidy fit for this " +
                         "sample, which is degenerate between ploidy basins. Its copy " +
                         "number is not comparable to samples with purity pinned."
            }

            // --- sex ----------------------------------------------------
            // Known from the clinical record. Supplying it stops mmkaryo
            // guessing from X and Y depth, and lets it flag a sample whose
            // sequenced sex does not match what the record says.
            def sex = (row.sex ?: '').trim().toUpperCase() ?: null
            if (sex && !(sex in ['XX', 'XY'])) {
                error "sample ${sample_id}: sex '${sex}' is not XX or XY. " +
                      "Leave the cell blank to infer it from depth."
            }

            // --- panel version ------------------------------------------
            // Each sample is analysed against the panel it was sequenced with.
            // Masking a v4 sample with the v7 panel discards genuine off-target
            // signal and leaves v4-only enriched regions unmasked, which
            // corrupts the copy-number profile in both directions.
            def panel_version = (row.panel_version ?: '').trim().toLowerCase() ?: null
            if (panel_version && !panelBedNames.containsKey(panel_version)) {
                error "sample ${sample_id}: unknown panel_version '${panel_version}'. " +
                      "Expected one of ${panelBedNames.keySet().join(', ')}."
            }

            def bed_t2t
            def bed_hg38
            if (panel_version) {
                bed_t2t  = file("${projectDir}/assets/${panelBedNames[panel_version].t2t}",
                                checkIfExists: true)
                bed_hg38 = file("${projectDir}/assets/${panelBedNames[panel_version].hg38}",
                                checkIfExists: true)
            }
            else {
                log.warn "sample ${sample_id}: no panel_version; falling back to the " +
                         "pipeline defaults. Confirm this is the panel the sample was " +
                         "actually sequenced with."
                bed_t2t  = file(params.panel_bed_t2t,  checkIfExists: true)
                bed_hg38 = file(params.panel_bed_hg38, checkIfExists: true)
            }

            def meta = [
                id            : sample_id,
                purity        : purity,
                purity_source : purity ? 'flow' : 'free_fit',
                sex           : sex,
                panel_version : panel_version,
                timepoint     : (row.timepoint ?: '').trim() ?: null,
                notes         : (row.notes ?: '').trim() ?: null,
            ]

            // checkIfExists: true so a bad path (file OR directory) fails fast
            // and legibly at launch, rather than deep inside MERGE_MINKNOW.
            tuple(meta,
                  file(row.minknow_bam.trim(), checkIfExists: true),
                  bed_t2t,
                  bed_hg38)
        }

    minknow_bams = ch_rows.map { meta, input, bed_t2t, bed_hg38 -> tuple(meta, input) }
    panel_beds   = ch_rows.map { meta, input, bed_t2t, bed_hg38 -> tuple(meta, bed_t2t, bed_hg38) }

    emit:
    minknow_bams = minknow_bams
    panel_beds   = panel_beds
}
