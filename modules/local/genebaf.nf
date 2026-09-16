process GENEBAF {
    tag      "${meta.id}"
    label    'process_single'

    input:
    tuple val(meta), path(bam), path(bai), path(panel_bed)

    output:
    tuple val(meta), path("${meta.id}.genebaf.tsv"),       emit: genes
    tuple val(meta), path("${meta.id}.genebaf.json"),      emit: summary
    tuple val(meta), path("${meta.id}.armloh.tsv"),        optional: true, emit: armloh
    tuple val(meta), path("${meta.id}.genebaf.sites.tsv"), emit: sites
    tuple val(meta), path("${meta.id}.genebaf.bins.tsv"),  emit: bins
    path "versions.yml",                                    emit: versions

    script:
    // Allelic state is read from on-target reads only.
    //
    // Off-target coverage here is 0.6 to 3x. At that depth a heterozygous site
    // shows one allele by chance about half the time, and a genotype-marginalised
    // likelihood has almost no power to separate balanced from imbalanced: 1 Mb
    // windows on a chromosome the depth caller reports as diploid returned minor
    // allele fractions anywhere from 0.01 to 0.50 with no centre, and widening to
    // 5 Mb did not help. Off-target BAF can see a complete loss of heterozygosity
    // and nothing else.
    //
    // On-target coverage is 17 to 57x. At that depth the unfolded allele fraction
    // at heterozygous sites reads the allelic ratio directly — on one sample the
    // panel windows on 1p were bimodal at 0.3 and 0.7, reproducing the depth
    // caller's 1p loss from independent reads.
    //
    // Loss of heterozygosity is tested per chromosome arm, not per gene. The
    // common SNPs inside one gene body sit in a handful of linkage blocks, so a
    // depleted heterozygous count there is far more often homozygosity across a
    // common haplotype than a somatic event: an earlier per-gene version put BRAF
    // in 14 of 21 samples, with 150 of 150 sites heterozygous in some and 0 of 150
    // in others, all at MAPQ 60. Pooling the windows on an arm gives observations
    // that do sit in different blocks.
    // genebaf_scope_span_v1: allelic state at the eight wide windows only;
    // depth on all. IGH/IGK/IGL are excluded inside genebaf.py regardless.
    def baf_regions_arg = params.genebaf_baf_regions ? "--baf-regions ${params.genebaf_baf_regions}" : ''
    """
    set -euo pipefail

    genebaf.py \\
        ${bam} \\
        ${panel_bed} \\
        ${params.baf_sites_vcf_hg38} \\
        ${meta.id} \\
        --sample ${meta.id} \\
        --cytobands ${params.cytoband_txt_hg38} \\
        --af-field ${params.genebaf_af_field} \\
        --min-depth ${params.genebaf_min_depth} \\
        --min-alt ${params.genebaf_min_alt} \\
        --min-mapq ${params.genebaf_min_mapq} \\
        --min-arm-windows ${params.genebaf_min_arm_windows} \\
        --min-arm-sites ${params.genebaf_min_arm_sites} \\
        --max-arm-p ${params.genebaf_max_arm_p} \\
        ${baf_regions_arg}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genebaf: \$(python3 -c "print('1.0')")
        pysam: \$(python3 -c "import pysam; print(pysam.__version__)")
    END_VERSIONS
    """

    stub:
    """
    printf 'gene\\tchrom\\tstart\\tend\\tmean_depth\\tn_sites\\tn_usable\\tn_het\\tn_hom\\tn_lowdepth\\tn_other_allele\\texpected_het_fraction\\tobserved_het_fraction\\tmedian_dev_from_half\\tfraction_in_central_band\\tallelic_state\\tstate_statistic\\n' > ${meta.id}.genebaf.tsv
    printf 'gene\\tchrom\\tpos\\tref_reads\\talt_reads\\talt_fraction\\n' > ${meta.id}.genebaf.sites.tsv
    printf 'gene\\tchrom\\tstart\\tend\\tmean_depth\\tlog2_ratio\\n' > ${meta.id}.genebaf.bins.tsv
    printf 'arm\\tchrom\\tn_windows\\tn_usable\\tn_het\\texpected_het_fraction\\tobserved_het_fraction\\tratio\\tp_value\\tverdict\\tgenes\\n' > ${meta.id}.armloh.tsv
    echo '{"sample": "${meta.id}", "states": {}, "arms": {}, "loh_arms": []}' > ${meta.id}.genebaf.json
    echo '"${task.process}": stub' > versions.yml
    """
}
