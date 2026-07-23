process BAF_LOH_SCREEN {
    tag      "cohort"
    label    'process_medium'

    input:
    // Cohort-scoped, unlike the per-sample processes elsewhere in the pipeline.
    // The screen normalises heterozygous site density for each region against
    // the cohort median for that same region, which cancels the confounding
    // effect of window size and mappability. That baseline needs at least three
    // samples, so the process consumes collected channels rather than running
    // once per sample.
    val  sample_ids               // list of sample identifiers, order matches the paths below
    path clair3_dirs, stageAs: 'clair3_*'
    path ichor_dirs,  stageAs: 'ichor_*'
    path panel_bed

    output:
    path "cohort.baf_screen.tsv",  emit: screen
    path "sample_map.tsv",         emit: sample_map
    path "ichor_map.tsv",          optional: true, emit: ichor_map
    path "versions.yml",           emit: versions

    script:
    def ids            = sample_ids instanceof List ? sample_ids : [sample_ids]
    def id_list        = ids.collect { it instanceof Map ? it.id : it }.join(' ')
    def min_site_depth = params.baf_loh_min_site_depth
    def min_sites      = params.baf_loh_min_sites
    def cdr            = params.baf_loh_cdr_threshold
    def bimodality     = params.baf_loh_bimodality_threshold
    def gap_distance   = params.baf_loh_max_gap_distance
    def use_ichor      = params.skip_ichorcna ? false : true
    """
    set -euo pipefail

    # Sample identifiers are written one per line rather than interpolated into
    # a loop, so that no shell metacharacter in an identifier can alter control
    # flow, and so the map is reproducible from the work directory alone.
    cat > sample_ids.txt <<'END_IDS'
${ids.collect { it instanceof Map ? it.id : it }.join('\n')}
END_IDS

    # Locate each sample's phased Clair3 VCF directory among the staged inputs.
    # Nextflow stages the collected directories as clair3_1, clair3_2, ...; the
    # sample each belongs to is identified by the phased VCF filenames inside
    # rather than by stage order, which is not guaranteed to match.
    : > sample_map.tsv
    while read -r sid; do
        [ -n "\$sid" ] || continue
        found=""
        for d in clair3_*; do
            [ -d "\$d" ] || continue
            # The phased per-chromosome VCFs live under the Clair3 tmp tree.
            cand="\$d/tmp/phase_output/phase_vcf"
            if [ -d "\$cand" ] && ls "\$cand"/phased_*.vcf.gz >/dev/null 2>&1; then
                # Confirm ownership by checking the sample column of any VCF.
                s_in_vcf=\$(bcftools query -l "\$(ls "\$cand"/phased_*.vcf.gz | head -1)" 2>/dev/null | head -1 || true)
                if [ "\$s_in_vcf" = "\$sid" ]; then
                    found="\$cand"
                    break
                fi
            fi
        done
        if [ -z "\$found" ]; then
            echo "WARNING: no phased VCF directory located for \$sid" >&2
            continue
        fi
        printf '%s\\t%s\\n' "\$sid" "\$found" >> sample_map.tsv
    done < sample_ids.txt

    if [ ! -s sample_map.tsv ]; then
        echo "ERROR: no samples resolved to a phased VCF directory" >&2
        ls -la >&2
        exit 1
    fi

    ICHOR_ARG=""
    if ${use_ichor}; then
        # Copy number files are passed as an explicit map. Staged directories do
        # not preserve the <dir>/<sample>/ichorcna_out/ layout the directory
        # based option expects, and matching on filename is unambiguous because
        # ichorCNA names its outputs after the sample.
        : > ichor_map.tsv
        while read -r sid _; do
            seg=\$(find -L ichor_* -name "\${sid}.cna.seg" -print -quit 2>/dev/null || true)
            par=\$(find -L ichor_* -name "\${sid}.params.txt" -print -quit 2>/dev/null || true)
            if [ -n "\$seg" ]; then
                printf '%s\\t%s\\t%s\\n' "\$sid" "\$seg" "\$par" >> ichor_map.tsv
            else
                echo "NOTE: no copy number segments for \$sid; BAF reported without a copy number call" >&2
            fi
        done < sample_map.tsv
        if [ -s ichor_map.tsv ]; then
            ICHOR_ARG="--ichor-map ichor_map.tsv"
        else
            # Copy number directories were supplied but nothing matched. This
            # is almost always a staging problem rather than genuinely absent
            # data, and it silently degrades every call to NO_CN, so it is
            # reported loudly rather than passed over.
            echo "WARNING: copy number directories were staged but no segment" >&2
            echo "         files matched any sample; screening BAF only." >&2
            ls -Ll ichor_* 2>/dev/null | head -20 >&2 || true
        fi
    fi

    baf_loh_screen.py \\
        --bed ${panel_bed} \\
        --sample-map sample_map.tsv \\
        \$ICHOR_ARG \\
        --min-site-depth ${min_site_depth} \\
        --min-sites ${min_sites} \\
        --cdr-threshold ${cdr} \\
        --bimodality-threshold ${bimodality} \\
        --max-gap-distance ${gap_distance} \\
        --out cohort.baf_screen.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        baf_loh_screen: v1.0
    END_VERSIONS
    """

    stub:
    """
    printf 'sample\\tregion\\tchrom\\tstart\\tend\\tflag\\n' > cohort.baf_screen.tsv
    touch sample_map.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
