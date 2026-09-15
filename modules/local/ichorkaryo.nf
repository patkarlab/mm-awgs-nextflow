process ICHORKARYO {
    tag      "${meta.id}"
    label    'process_single'

    input:
    tuple val(meta), path(seg), path(params_file), path(cna_seg), path(wig), path(panel_bed)

    output:
    tuple val(meta), path("${meta.id}.ichorkaryo.json"), emit: karyotype
    path "versions.yml",                                 emit: versions

    script:
    // Depth-based copy number comes from ichorCNA rather than mmkaryo.
    //
    // mmkaryo masks the adaptive sampling panel before segmenting, and on this
    // cohort the panel window over TP53 covers the region a del(17p) FISH probe
    // targets — so a focal 17p deletion is removed from the data before the
    // caller sees it. On the sample with FISH monosomy 13 and del(17p) both at
    // 95 percent, mmkaryo returned a near-flat genome; ichorCNA calls chr13 at
    // copy number 1 over 93 bins and chr17p over 6, because it bins the whole
    // genome at 1 Mb without masking.
    //
    // Purity is supplied from flow cytometry, not taken from ichorCNA. What
    // ichorCNA reports as tumour fraction is the magnitude of copy-number
    // deviation, which is the right quantity for cfDNA and the wrong one for
    // sorted plasma cells: four samples here returned a fraction near zero
    // while carrying driver variants at 35 to 69 percent allele fraction. That
    // value travels through as cna_burden, named for what it measures.
    //
    // --cna and --wig feed the qc block. The per-bin log2 in the .cna.seg is
    // what the detection limits are computed from: the segment-level sigma
    // already in use for the event threshold is the scatter of segment
    // medians, which is far smaller and answers a different question. The WIG
    // supplies read counts, which ichorCNA's own outputs do not carry.
    //
    // Both are staged as an empty list when the upstream channel had nothing
    // for this sample, which is falsy, so the flag is simply omitted and the
    // corresponding qc fields come back null.
    def purity_arg = meta.purity ? "--purity ${meta.purity}" : ''
    def sex_arg    = meta.sex    ? "--sex ${meta.sex}"       : ''
    def params_arg = params_file ? "--params ${params_file}" : ''
    def cna_arg    = cna_seg     ? "--cna ${cna_seg}"        : ''
    def wig_arg    = wig         ? "--wig ${wig}"            : ''
    """
    set -euo pipefail

    ichorkaryo.py \\
        --seg ${seg} \\
        ${params_arg} \\
        ${cna_arg} \\
        ${wig_arg} \\
        --cytobands ${params.cytoband_txt_hg38} \\
        --panel-bed ${panel_bed} \\
        --sample ${meta.id} \\
        ${purity_arg} \\
        ${sex_arg} \\
        --min-log2 ${params.ichorkaryo_min_log2} \\
        --n-sigma ${params.ichorkaryo_n_sigma} \\
        --out ${meta.id}.ichorkaryo.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ichorkaryo: \$(python3 -c "print('1.0')")
        python: \$(python3 --version | sed 's/Python //')
    END_VERSIONS
    """

    stub:
    """
    echo '{"sample": "${meta.id}", "arms": {}, "genes": [], "altered_segments": []}' \\
        > ${meta.id}.ichorkaryo.json
    echo '"${task.process}": stub' > versions.yml
    """
}
