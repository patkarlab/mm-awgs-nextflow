process MMPLOT {
    tag      "${meta.id}"
    label    'process_low'

    input:
    tuple val(meta), path(bins), path(segments_baf), path(windows), path(karyotype_json), path(panel_bed)
    path cytobands
    path genes

    output:
    tuple val(meta), path("${meta.id}.cn_genome.png"),  emit: genome
    tuple val(meta), path("${meta.id}.cn_grid.png"),    emit: grid
    tuple val(meta), path("${meta.id}.cn_chr*.png"),    optional: true, emit: chromosomes
    path "versions.yml",                                emit: versions

    script:
    def cytoband_arg = cytobands.name != 'NO_FILE' ? "--cytobands ${cytobands}" : ''
    def gene_arg     = genes.name     != 'NO_FILE' ? "--genes ${genes}"         : ''
    def target_arg   = panel_bed.name != 'NO_FILE' ? "--targets ${panel_bed}"   : ''
    """
    set -euo pipefail

    # Copy number cannot be computed without a purity, so take whichever the
    # upstream fit settled on rather than re-deriving it here. The subtitle
    # carries the ISCN string so the figure is self-describing.
    PURITY=\$(python -c "import json; print(json.load(open('${karyotype_json}'))['fit']['purity'])")
    PLOIDY=\$(python -c "import json; print(json.load(open('${karyotype_json}'))['fit']['ploidy'])")
    ISCN=\$(python -c "import json; print(json.load(open('${karyotype_json}')).get('ISCN_karyotype',''))")

    COMMON=(--purity "\$PURITY" --ploidy "\$PLOIDY" --windows ${windows}
            --max-cn ${params.mmkaryo_max_cn} --title "${meta.id}"
            ${cytoband_arg} ${gene_arg} ${target_arg})

    # Genome overview and the compact grid.
    mmplot.py ${bins} ${segments_baf} ${meta.id}.cn_genome.png \\
        "\${COMMON[@]}" --subtitle "\$ISCN" --genome
    mmplot.py ${bins} ${segments_baf} ${meta.id}.cn_grid.png \\
        "\${COMMON[@]}" --subtitle "\$ISCN"

    # One page per chromosome. The grid is a genome overview; these are what a
    # reader opens to look at a specific event, with gene highlight bands, the
    # panel targets track and a labelled cytogram at full width.
    # Chromosomes absent from the segment table are skipped rather than
    # producing an empty page, which is why mmplot exits non-zero there.
    for C in ${params.mmplot_chromosomes}; do
        mmplot.py ${bins} ${segments_baf} ${meta.id}.cn_chr\${C}.png \\
            "\${COMMON[@]}" --chr chr\${C} || \\
            echo "  no data for chr\${C}, skipped" >&2
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cn_genome.png ${meta.id}.cn_grid.png versions.yml
    """
}
