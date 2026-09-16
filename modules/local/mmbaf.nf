// RETIRED 2026-09-16 (retire_mmkaryo_v1). No longer included by any
// workflow. Copy-number figures come from ICHORCNA_PLOT on the hg38 track;
// see modules/local/ichorcna_plot.nf. Kept for one release so the history
// of the copy-number tab stays readable, then delete.
process MMBAF {
    tag      "${meta.id}"
    label    'process_high'

    input:
    tuple val(meta), path(bam), path(bai), path(segments), path(karyotype_json)
    path t2t_fai
    path g1000_vcf

    output:
    tuple val(meta), path("${meta.id}.segments.baf.tsv"), emit: segments
    tuple val(meta), path("${meta.id}.baf_windows.tsv"),  emit: windows
    tuple val(meta), path("${meta.id}.baf.json"),         emit: summary
    path "versions.yml",                                  emit: versions

    script:
    // The joint purity/ploidy refit is only well posed with purity pinned.
    // Left free it can reinterpret the whole genome, because (3,1) at purity
    // 0.95 and (5,1) at purity 0.33 produce nearly the same allele fraction.
    // Prefer the sample sheet value; fall back to mmkaryo's depth fit.
    def purity_arg = meta.purity ? "--purity ${meta.purity}" : ''
    """
    set -euo pipefail

    EXTRA="${purity_arg}"
    if [[ -z "\$EXTRA" ]]; then
        DEPTH_PURITY=\$(python -c "import json; print(json.load(open('${karyotype_json}'))['fit']['purity'])" 2>/dev/null || true)
        if [[ -n "\$DEPTH_PURITY" && "\$DEPTH_PURITY" != "None" ]]; then
            EXTRA="--depth-purity \$DEPTH_PURITY"
        fi
    fi

    # --af-field selects the population allele frequency used as the
    # heterozygosity prior. Match it to the cohort; SAS_AF is correct here and
    # measurably sharper than the global AF.
    #
    # --jobs controls wall clock: pysam.count_coverage is effectively the whole
    # runtime and chromosomes are independent, so this is ~12x on 24 workers.
    mmbaf.py \\
        ${bam} \\
        ${t2t_fai} \\
        ${g1000_vcf} \\
        ${segments} \\
        ${meta.id} \\
        --af-field ${params.mmbaf_af_field} \\
        --min-depth ${params.mmbaf_min_depth} \\
        --min-sites ${params.mmbaf_min_sites} \\
        --window ${params.mmbaf_window} \\
        --jobs ${task.cpus} \\
        --threads 2 \\
        \$EXTRA \\
        > ${meta.id}.baf.stdout.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //')
        numpy: \$(python -c "import numpy; print(numpy.__version__)")
        pysam: \$(python -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.segments.baf.tsv ${meta.id}.baf_windows.tsv
    echo '{}' > ${meta.id}.baf.json
    touch versions.yml
    """
}
