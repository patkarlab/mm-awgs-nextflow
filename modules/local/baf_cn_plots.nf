process BAF_CN_PLOTS {
    tag      "cohort"
    label    'process_low'

    input:
    path screen_tsv
    path sample_map
    path clair3_dirs, stageAs: 'clair3_*'
    path ichor_dirs,  stageAs: 'ichor_*'
    path panel_bed

    output:
    path "baf_cn_figures", emit: figures
    path "versions.yml",   emit: versions

    script:
    def min_site_depth = params.baf_loh_min_site_depth
    def min_sites      = params.baf_loh_min_sites
    """
    set -euo pipefail
    mkdir -p baf_cn_figures

    # The plotting script takes a four-column table. It is rebuilt here from the
    # sample map the screen already resolved, so the figures are guaranteed to
    # be drawn from exactly the inputs that were screened.
    printf 'sample\\tvcf_path\\tcna_seg\\tparams_txt\\n' > plot_inputs.tsv
    while read -r sid vcfdir; do
        [ -n "\$sid" ] || continue
        seg=\$(find -L ichor_* -name "\${sid}.cna.seg" -print -quit 2>/dev/null || true)
        par=\$(find -L ichor_* -name "\${sid}.params.txt" -print -quit 2>/dev/null || true)
        printf '%s\\t%s\\t%s\\t%s\\n' "\$sid" "\$vcfdir" "\${seg:-NA}" "\${par:-NA}" >> plot_inputs.tsv
    done < ${sample_map}

    plot_baf_cn.py \\
        --screen ${screen_tsv} \\
        --bed ${panel_bed} \\
        --sample-map plot_inputs.tsv \\
        --min-site-depth ${min_site_depth} \\
        --min-sites ${min_sites} \\
        --outdir baf_cn_figures

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p baf_cn_figures
    touch baf_cn_figures/.stub
    echo '"${task.process}": stub' > versions.yml
    """
}
