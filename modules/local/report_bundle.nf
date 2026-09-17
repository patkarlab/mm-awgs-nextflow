process REPORT_BUNDLE {
    tag      "bundle"
    label    'process_low'

    input:
    // These inputs are consumed purely to establish ordering. The bundle is
    // assembled by scanning the published output tree, which only exists once
    // publishDir has run for every producing process, so the process must not
    // start until those have completed. Nextflow provides no direct "after
    // publish" dependency, so the completion signals are taken as inputs and
    // ignored in the script body.
    val  ready_signals
    // bundle_guard_v1: every sample in the run, and the artefact kinds the
    // current configuration must have published for each of them.
    val  sample_ids
    val  required_kinds
    path bundle_script
    // The panel BEDs are staged rather than resolved relative to the
    // script, because inside a work directory bin/ is staged alone with no
    // assets/ beside it. Without them the alignment slicing skips silently
    // and the bundle ships without the reads its calls were made from.
    path panel_bed_t2t
    path panel_bed_hg38
    // allelic_gene_track_v1: hg38 gene model, copied into each sample's
    // allelic/ so the window plots can draw gene boundaries. May be an
    // empty list when params.gene_model_hg38 is null.
    path gene_model_hg38
    val  outdir
    val  bundle_name

    output:
    path "${bundle_name}",        emit: bundle
    path "${bundle_name}.zip",    optional: true, emit: archive
    path "${bundle_name}.tar.gz", optional: true, emit: tarball
    path "versions.yml",          emit: versions

    script:
    // The results directory is resolved to an absolute path because the script
    // runs from the task work directory, not the launch directory.
    def results_dir = outdir.toString().startsWith('/')
        ? outdir.toString()
        : "${workflow.launchDir}/${outdir}"
    """
    set -euo pipefail

    export PANEL_BED_T2T="\$(readlink -f ${panel_bed_t2t})"
    export PANEL_BED_HG38="\$(readlink -f ${panel_bed_hg38})"
    export GENE_MODEL_HG38="${gene_model_hg38 ? '\$(readlink -f ' + gene_model_hg38 + ')' : ''}"
    # Sliced BAMs are excluded by default. A panel slice is not a reduction
    # under adaptive sampling, where nearly all depth is on-target: the
    # 21-sample cohort bundle reached 13 GB, and zipping BGZF gains nothing.
    # The IGV pages already carry embedded reads at each call.
    export BUNDLE_NO_BAMS="${params.bundle_no_bams ? '1' : ''}"

    if [ ! -d "${results_dir}" ]; then
        echo "ERROR: results directory not found: ${results_dir}" >&2
        echo "The bundle is assembled from published outputs, so it cannot run" >&2
        echo "before publishDir has written them." >&2
        exit 1
    fi

    # bundle_guard_v1. Completion of the producing tasks is already
    # guaranteed by ready_signals. What is not guaranteed is that publishDir
    # has finished copying their outputs into the results tree, and on a
    # resumed run every cached task re-publishes at once. Poll until the
    # tree holds every required artefact for every sample, then fail
    # loudly rather than ship a bundle with silent gaps.
    printf '%s\\n' ${sample_ids.join(' ')} > sample_ids.txt
    deadline=\$(( SECONDS + ${params.bundle_publish_wait_s} ))
    until check_published_outputs.sh "${results_dir}" sample_ids.txt "${required_kinds}" \\
            > publish_check.log 2>&1; do
        if (( SECONDS >= deadline )); then
            cat publish_check.log >&2
            echo "ERROR: required outputs still absent from ${results_dir} after" >&2
            echo "       ${params.bundle_publish_wait_s} s. Either publishDir has not" >&2
            echo "       written them or a producing process did not emit them." >&2
            echo "       Re-run bin/check_published_outputs.sh by hand to see which." >&2
            exit 1
        fi
        sleep 30
    done
    cat publish_check.log

    export BUNDLE_STRICT="${params.bundle_strict ? '1' : '0'}"
    export BUNDLE_REQUIRE="${required_kinds}"

    bash ${bundle_script} "${results_dir}" "${bundle_name}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        build_report_bundle: v1.0
        zip: \$( { zip --version 2>/dev/null || true; } | sed -n '2s/.*Zip \\([0-9.]*\\).*/\\1/p' || echo "not available")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${bundle_name}
    touch ${bundle_name}/.stub
    echo '"${task.process}": stub' > versions.yml
    """
}
