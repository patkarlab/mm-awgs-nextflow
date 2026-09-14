process NANOMONSV {
    tag      "${meta.id}"
    label    'process_medium'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.id}.nanomonsv.result.txt"), emit: result
    tuple val(meta), path("${meta.id}.nanomonsv.result.vcf"), emit: vcf
    tuple val(meta), path("${meta.id}.nanomonsv.sbnd.result.txt"), optional: true, emit: sbnd
    path "versions.yml",                                      emit: versions

    script:
    // Tumour-only, with a non-matched control panel doing the germline and
    // common-SV filtering. There are no matched normals for this cohort.
    //
    // Why this runs alongside the four-caller ensemble rather than inside it
    // ---------------------------------------------------------------------
    // nanomonsv is far more specific than the ensemble and less sensitive at
    // the support levels this cohort operates at. Measured on three samples:
    //
    //   A cytogenetically confirmed IGH translocation at moderate support was
    //   recovered with both reciprocal breakends, PASS, from 44 PASS calls
    //   genome-wide against roughly 950 rows for the same sample from the
    //   merged ensemble.
    //
    //   On a sample with no structural findings, three interchromosomal calls
    //   survived, all at the minimum supporting-read threshold, and none in an
    //   Ig locus. The ensemble reported dozens of Ig-partnered junctions for
    //   the same sample.
    //
    //   A second confirmed translocation at two supporting reads was captured
    //   by the parse step at MAPQ 60 on both breakends but never reported.
    //   Relaxing min_tumor_variant_read_num, min_tumor_VAF,
    //   max_overhang_size_thres and cluster_margin_size did not recover it, so
    //   the loss is in cluster formation rather than any exposed threshold.
    //
    // The control panel is what produces the specificity: V(D)J recombination
    // and somatic hypermutation are present in every normal B cell, and a
    // panel of normals subtracts them. No trained classifier in the stack
    // manages that.
    //
    // So its result is a confirmation layer joined onto the annotated table,
    // not a fifth vote in the SURVIVOR merge. A junction the ensemble reports
    // and nanomonsv confirms is strongly supported; one nanomonsv does not
    // report means nothing either way, because two-read junctions are below
    // what it will emit.
    //
    // --qv20 matches dorado basecalling on R10.4.1 chemistry. The shipped
    // thresholds are lower than the tool's defaults (3 reads, 0.05 VAF)
    // because adaptive sampling puts real translocation breakpoints at
    // single-digit support; they are parameters so the choice stays visible.
    def repeat_arg = params.nanomonsv_simple_repeat_bed ?
        "--simple_repeat_bed ${params.nanomonsv_simple_repeat_bed}" : ''
    def panel_arg = params.nanomonsv_control_panel ?
        "--control_panel_prefix ${params.nanomonsv_control_panel}" : ''
    """
    set -euo pipefail

    nanomonsv parse \\
        ${bam} \\
        ${meta.id}

    nanomonsv get \\
        ${meta.id} \\
        ${bam} \\
        ${params.t2t_fasta} \\
        ${panel_arg} \\
        ${repeat_arg} \\
        --qv20 \\
        --processes ${task.cpus} \\
        --min_tumor_variant_read_num ${params.nanomonsv_min_read_num} \\
        --min_tumor_VAF ${params.nanomonsv_min_vaf}

    if [ ! -s "${meta.id}.nanomonsv.result.txt" ]; then
        echo "ERROR: nanomonsv produced no result file for ${meta.id}" >&2
        ls -la >&2 || true
        exit 1
    fi

    n_pass=\$(awk -F'\\t' 'NR>1 && \$NF=="PASS"' ${meta.id}.nanomonsv.result.txt | wc -l)
    n_inter=\$(awk -F'\\t' 'NR>1 && \$NF=="PASS" && \$1!=\$4' ${meta.id}.nanomonsv.result.txt | wc -l)
    echo "nanomonsv ${meta.id}: \$n_pass PASS, \$n_inter interchromosomal"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanomonsv: \$(nanomonsv --version 2>&1 | grep -o '[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' | head -1)
    END_VERSIONS
    """

    stub:
    """
    printf 'Chr_1\\tPos_1\\tDir_1\\tChr_2\\tPos_2\\tDir_2\\tInserted_Seq\\tSV_ID\\tChecked_Read_Num_Tumor\\tSupporting_Read_Num_Tumor\\tSupporting_Read_Num_Tumor_HP_BP1\\tSupporting_Read_Num_Tumor_HP_BP2\\tIs_Filter\\n' \\
        > ${meta.id}.nanomonsv.result.txt
    touch ${meta.id}.nanomonsv.result.vcf
    echo '"${task.process}": stub' > versions.yml
    """
}
