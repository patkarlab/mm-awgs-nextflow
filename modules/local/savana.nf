process SAVANA {
    tag      "${meta.id}"
    label    'process_high'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    // classified.vcf, NOT classified.somatic.vcf. See the note below: the
    // somatic model discards the calls this project exists to detect.
    tuple val(meta), path("${meta.id}.savana.vcf"),          emit: vcf
    tuple val(meta), path("${meta.id}.savana.somatic.vcf"),  emit: somatic_vcf
    tuple val(meta), path("savana_out"),                      emit: outdir
    path "versions.yml",                                       emit: versions

    script:
    // Why the unfiltered classified.vcf is the caller output
    // ------------------------------------------------------
    // SAVANA's ont-somatic-to.pkl tumour-only model was trained on
    // whole-genome data, where a junction with two supporting reads against
    // 50x depth is noise. Under adaptive sampling at a breakpoint cluster
    // region, two reads can be the entire depth available, and the model has
    // no basis for that feature vector.
    //
    // Measured on this cohort's cytogenetically confirmed IGH translocation
    // cases: each was recovered by the clustering step with both breakends at
    // MAPQ 60, matching an independent caller's coordinates to within one
    // base, and each was then classified PREDICTED_NOISE. On one such sample
    // no call the model passed sat below six supporting reads, while the
    // majority of its rejected calls sat at two. The model behaves as a depth
    // filter on this data, which is the same floor that caused a
    // default-parameter caller to miss those junctions in the first place.
    //
    // So SAVANA is used for its clustering, and its classification travels as
    // an annotation. CLASS reaches the annotated TSV through the VCF FILTER
    // column so a reviewer can see the verdict; it does not gate anything.
    //
    // The somatic VCF is still emitted, unused by the ensemble, so the two can
    // be compared directly.
    def snp_arg = params.savana_snp_resource ? "-vg ${params.savana_snp_resource}" : ''
    """
    set -euo pipefail
    mkdir -p savana_out

    savana to \\
        --tumour ${bam} \\
        --ref ${params.t2t_fasta} \\
        --outdir savana_out \\
        --sample ${meta.id} \\
        --ont \\
        --threads ${task.cpus} \\
        --min_support ${params.savana_min_support} \\
        --overwrite \\
        ${snp_arg}

    # Normalize both outputs to predictable names. SAVANA prefixes with the
    # --sample value, which is meta.id here, but the suffix has varied between
    # versions so each candidate is tried rather than assumed.
    for pair in \\
        "${meta.id}.savana.vcf:savana_out/${meta.id}.classified.vcf" \\
        "${meta.id}.savana.somatic.vcf:savana_out/${meta.id}.classified.somatic.vcf"
    do
        dest="\${pair%%:*}"
        src="\${pair#*:}"
        if [ -s "\$src" ]; then
            cp "\$src" "\$dest"
        else
            # Fall back to a glob in case the sample prefix differs.
            alt=\$(find savana_out -maxdepth 1 -name "*\$(basename "\$src" | sed 's/^[^.]*//')" | head -1)
            if [ -n "\$alt" ] && [ -s "\$alt" ]; then
                cp "\$alt" "\$dest"
            else
                echo "ERROR: SAVANA produced no \$(basename "\$src") for ${meta.id}" >&2
                ls -lR savana_out/ >&2 || true
                exit 1
            fi
        fi
    done

    n_all=\$(grep -vc '^#' ${meta.id}.savana.vcf || true)
    n_som=\$(grep -vc '^#' ${meta.id}.savana.somatic.vcf || true)
    echo "SAVANA ${meta.id}: \$n_all classified calls, \$n_som passing the somatic model"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        savana: \$(savana --version 2>&1 | grep -o '[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    mkdir -p savana_out
    touch ${meta.id}.savana.vcf ${meta.id}.savana.somatic.vcf
    echo '"${task.process}": stub' > versions.yml
    """
}
