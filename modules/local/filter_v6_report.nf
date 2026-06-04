process FILTER_V6_REPORT {
    tag      "${meta.id}"
    label    'process_low'

    input:
    // From VEP_ANNOTATE_CLAIR3.out.candidates_tsv
    tuple val(meta), path(candidates_tsv)

    output:
    tuple val(meta), path("v6_filtered/*.v6_clinical.tsv"), optional: true, emit: clinical
    tuple val(meta), path("v6_filtered/*.v6_filtered.tsv"), optional: true, emit: filtered
    tuple val(meta), path("v6_filtered/v6_filter_summary.tsv"), optional: true, emit: summary
    path "versions.yml", emit: versions

    script:
    // Gene-symbol panel filter + protein-altering consequence filter + exon
    // date-fix + AF% + count passthrough. Pure stdlib python (no pandas), so it
    // runs in the awgs_sv env. REF_COUNT/ALT_COUNT/DP are emitted by
    // VEP_ANNOTATE_CLAIR3 (from Clair3 AD/DP) and passed through here.
    def include_ig = params.v6_include_ig ? '--include-ig' : ''
    """
    filter_v6_somatic_candidates.py \\
        --input ${candidates_tsv} \\
        --outdir v6_filtered \\
        --max-pop-af ${params.v6_max_pop_af} \\
        ${include_ig}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_v6_somatic_candidates.py: v1-stdlib
    END_VERSIONS
    """

    stub:
    """
    mkdir -p v6_filtered
    touch v6_filtered/${meta.id}.v6_clinical.tsv
    touch v6_filtered/${meta.id}.v6_filtered.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
