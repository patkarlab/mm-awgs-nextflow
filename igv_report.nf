process IGV_REPORT {
    tag      "${meta.id}:${mode}"
    label    'process_medium'

    // Renders one self-contained igv-reports HTML for a single sample, for one
    // evidence class (mode = 'somatic' or 'translocations'). The same process
    // serves both tracks; the caller passes the matching sites TSV, BAM, and
    // reference. Uses bin/igv_snapshots.py (stdlib-only, awgs_sv).
    //
    // mode = 'somatic'        sites = v6 clinical TSV,      bam = hg38 BAM
    // mode = 'translocations' sites = translocations TSV,   bam = T2T  BAM

    input:
    // meta.id is the sequencing id only (no patient identifiers).
    tuple val(meta), val(mode), path(sites_tsv), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), val(mode), path("igv_out/${meta.id}.${mode}.html"), emit: html
    path "versions.yml",                                                 emit: versions

    script:
    def flanking = params.igv_flanking ?: 5000
    """
    set -euo pipefail
    mkdir -p igv_out

    # igv_snapshots.py handles defensive (re)indexing, the legitimate
    # zero-sites case (writes a placeholder HTML, exits 0), and surfaces the
    # full create_report command on failure.
    igv_snapshots.py \\
        --mode ${mode} \\
        --sample ${meta.id} \\
        --sites-tsv ${sites_tsv} \\
        --bam ${bam} \\
        --fasta ${fasta} \\
        --out-html igv_out/${meta.id}.${mode}.html \\
        --flanking ${flanking}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        igv_snapshots: stdlib
        igv_reports: \$(create_report --version 2>/dev/null || echo "unknown")
        samtools: \$(samtools --version | head -1 | awk '{print \$NF}')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p igv_out
    touch igv_out/${meta.id}.${mode}.html
    echo '"${task.process}": stub' > versions.yml
    """
}
