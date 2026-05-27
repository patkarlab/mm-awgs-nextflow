// Two named processes sharing the same logic. The split exists so each
// reference's index publishes under its own output subdirectory and shows
// up cleanly in the Nextflow DAG.

process SAMTOOLS_INDEX_T2T {
    tag      "${meta.id}"
    label    'process_low'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: bam_bai
    path "versions.yml",                             emit: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai
    echo '"${task.process}": stub' > versions.yml
    """
}

process SAMTOOLS_INDEX_HG38 {
    tag      "${meta.id}"
    label    'process_low'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path(bam), path("${bam}.bai"), emit: bam_bai
    path "versions.yml",                             emit: versions

    script:
    """
    samtools index -@ ${task.cpus} ${bam}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${bam}.bai
    echo '"${task.process}": stub' > versions.yml
    """
}
