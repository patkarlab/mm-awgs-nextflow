process SURVIVOR_MERGE {
    tag      "${meta.id}"
    label    'process_low'

    input:
    // Four caller VCFs in fixed order: Sniffles, CuteSV, Severus, SAVANA.
    // Sniffles and CuteSV are .vcf.gz; Severus and SAVANA are plain .vcf.
    // SAVANA is LAST so the first three SUPP_VEC bit positions are unchanged.
    tuple val(meta), path(sniffles_vcf_gz), path(cutesv_vcf_gz), path(severus_vcf), path(savana_vcf)

    output:
    tuple val(meta), path("${meta.id}.merged.vcf.gz"), path("${meta.id}.merged.vcf.gz.tbi"), emit: merged_vcf
    path "versions.yml",                                                                       emit: versions

    script:
    """
    set -euo pipefail

    # SURVIVOR consumes plain (uncompressed) VCFs only.
    zcat ${sniffles_vcf_gz} > sniffles.vcf
    zcat ${cutesv_vcf_gz}   > cutesv.vcf
    if [[ "${severus_vcf}" == *.gz ]]; then
        zcat ${severus_vcf} > severus.vcf
    else
        cp ${severus_vcf} severus.vcf
    fi
    if [[ "${savana_vcf}" == *.gz ]]; then
        zcat ${savana_vcf} > savana.vcf
    else
        cp ${savana_vcf} savana.vcf
    fi

    # The caller order written to vcflist.txt is the order SURVIVOR uses for
    # SUPP_VEC bits. We MUST keep this exact order — downstream MM annotation
    # decodes 1000/0100/0010/0001 as sniffles/cutesv/severus/savana
    # respectively. SAVANA was appended last precisely so the first three
    # positions still mean what they meant when the vector was three bits.
    cat > vcflist.txt <<'EOF'
sniffles.vcf
cutesv.vcf
severus.vcf
savana.vcf
EOF

    # SURVIVOR merge args (matches production survivor_merge_t2t.sh):
    #   max_dist_bp     ${params.survivor_max_dist}
    #   min_callers     ${params.survivor_min_callers}
    #   take_type       ${params.survivor_take_type}    same SV type required
    #   take_strand     ${params.survivor_take_strand}  strand match required
    #   estimate_dist   ${params.survivor_estimate_dist} 0 = use actual positions
    #   min_size        ${params.survivor_min_size}
    SURVIVOR merge \\
        vcflist.txt \\
        ${params.survivor_max_dist} \\
        ${params.survivor_min_callers} \\
        ${params.survivor_take_type} \\
        ${params.survivor_take_strand} \\
        ${params.survivor_estimate_dist} \\
        ${params.survivor_min_size} \\
        ${meta.id}.merged.vcf

    # SURVIVOR's output is sorted at chromosome granularity but not strictly
    # position-sorted within chromosomes. bcftools sort handles this AND
    # bgzips in one pass.
    mkdir -p bcftools_tmp
    bcftools sort \\
        -O z \\
        -o ${meta.id}.merged.vcf.gz \\
        --temp-dir bcftools_tmp \\
        ${meta.id}.merged.vcf

    tabix -p vcf -f ${meta.id}.merged.vcf.gz

    # Cleanup intermediate plain VCFs
    rm -f sniffles.vcf cutesv.vcf severus.vcf savana.vcf ${meta.id}.merged.vcf
    rm -rf bcftools_tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SURVIVOR: \$( { SURVIVOR --version 2>&1 || true; } | head -1 | awk '{print \$NF}' || echo unknown)
        bcftools: \$( { bcftools --version || true; } | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.merged.vcf.gz ${meta.id}.merged.vcf.gz.tbi
    echo '"${task.process}": stub' > versions.yml
    """
}
