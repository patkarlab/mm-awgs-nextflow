process ANNOTATE_MM_TRANSLOCATIONS {
    tag      "${meta.id}"
    label    'process_low'

    input:
    tuple val(meta), path(merged_vcf), path(merged_tbi)

    output:
    tuple val(meta), path("${meta.id}.mm_annotated.tsv"), emit: tsv
    path "versions.yml",                                    emit: versions

    script:
    // Each of the three is passed only when its parameter is set, so
    // the annotator degrades to naming breakpoints from the panel
    // interval rather than failing. Without the gene model the
    // distance columns stay empty and compound intervals keep their
    // compound labels; without the anchors only dictionary-named
    // pairs are reportable; without the exclusion list nothing is
    // dropped, which is the shipped default.
    def gene_model = params.gene_model_t2t
        ? "--gene-model ${params.gene_model_t2t}" : ''
    def ig_segments = params.ig_segments_t2t
        ? "--ig-segments ${params.ig_segments_t2t}" : ''
    def anchors = params.mm_translocation_anchors
        ? "--anchors ${params.mm_translocation_anchors}" : ''
    def excluded = params.mm_excluded_junctions
        ? "--excluded-junctions ${params.mm_excluded_junctions}" : ''
    """
    annotate_mm_translocations.py \\
        --vcf ${merged_vcf} \\
        --panel-bed ${params.panel_bed_t2t} \\
        --cytoband-bed ${params.cytoband_bed_t2t} \\
        --dictionary ${params.mm_translocation_dict} \\
        ${gene_model} \\
        ${anchors} \\
        ${ig_segments} \\
        ${excluded} \\
        --sample ${meta.id} \\
        --output ${meta.id}.mm_annotated.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        annotate_mm_translocations.py: \$(annotate_mm_translocations.py --version 2>&1 || echo 'v0.1')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.mm_annotated.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
