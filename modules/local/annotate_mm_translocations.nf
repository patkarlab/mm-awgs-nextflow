process ANNOTATE_MM_TRANSLOCATIONS {
    tag      "${meta.id}"
    label    'process_low'

    input:
    tuple val(meta), path(merged_vcf), path(merged_tbi)

    output:
    tuple val(meta), path("${meta.id}.mm_annotated.tsv"), emit: tsv
    path "versions.yml",                                    emit: versions

    script:
    """
    annotate_mm_translocations.py \\
        --vcf ${merged_vcf} \\
        --panel-bed ${params.panel_bed_t2t} \\
        --cytoband-bed ${params.cytoband_bed_t2t} \\
        --dictionary ${params.mm_translocation_dict} \\
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
