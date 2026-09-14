process MMKARYO {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    // panel_bed rides on the sample tuple: the cohort spans four panel
    // versions and each sample must be masked with the panel it was
    // sequenced with, not a pipeline-wide default.
    tuple val(meta), path(bam), path(bai), path(panel_bed)
    path t2t_fai
    path t2t_fasta
    path cytobands

    output:
    tuple val(meta), path("${meta.id}.bins.tsv"),         emit: bins
    tuple val(meta), path("${meta.id}.segments.tsv"),     emit: segments
    tuple val(meta), path("${meta.id}.arms.tsv"),         emit: arms
    tuple val(meta), path("${meta.id}.karyotype.json"),   emit: karyotype
    tuple val(meta), path("${meta.id}.cytobands.tsv"),    optional: true, emit: cytobands
    tuple val(meta), path("${meta.id}.karyotype.png"),    optional: true, emit: plot
    path "versions.yml",                                  emit: versions

    script:
    // Purity comes from the sample sheet. Depth alone cannot separate purity
    // from ploidy, so leaving it free lands the fit in a low-purity basin on
    // some samples: two cohort cases fitted 0.26 and 0.30 against measured
    // flow values of 0.88 and 0.98. Supply it wherever it is known.
    def purity_arg    = meta.purity    ? "--purity ${meta.purity}"        : ''
    // Sex comes from the sample sheet when known. Absent, mmkaryo infers it
    // from X and Y depth, which is what it did for the whole cohort before.
    def sex_arg       = meta.sex       ? "--sex ${meta.sex}"                : ''
    def cytoband_arg  = cytobands.name != 'NO_FILE' ? "--cytobands ${cytobands}" : ''
    def reference_arg = t2t_fasta.name != 'NO_FILE' ? "--reference ${t2t_fasta}" : ''
    """
    set -euo pipefail

    # The adaptive sampling panel is masked, not corrected for. Panel
    # boundaries would otherwise create copy-number changepoints, and masking
    # removes the need for a panel of normals in the genome-wide tier.
    # --mask-stains additionally drops acrocentric short arms, centromeres and
    # heterochromatic variants, where rDNA and satellite make read counts
    # meaningless: without it 13p, 14p, 15p, 21p and 22p return log2 -2 to -8,
    # which is "nothing maps" rather than "deleted".
    mmkaryo.py \\
        ${bam} \\
        ${t2t_fai} \\
        ${meta.id} \\
        --exclude ${panel_bed} \\
        --bin-size ${params.mmkaryo_bin_size} \\
        --min-mapq ${params.mmkaryo_min_mapq} \\
        --min-delta ${params.mmkaryo_min_delta} \\
        --min-seg-bins ${params.mmkaryo_min_seg_bins} \\
        --max-cn ${params.mmkaryo_max_cn} \\
        --threads ${task.cpus} \\
        --sample ${meta.id} \\
        ${reference_arg} \\
        ${cytoband_arg} \\
        ${purity_arg} \\
        ${sex_arg} \\
        > ${meta.id}.karyotype.stdout.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.bins.tsv ${meta.id}.segments.tsv ${meta.id}.arms.tsv
    touch ${meta.id}.cytobands.tsv ${meta.id}.karyotype.png
    echo '{}' > ${meta.id}.karyotype.json
    touch versions.yml
    """
}
