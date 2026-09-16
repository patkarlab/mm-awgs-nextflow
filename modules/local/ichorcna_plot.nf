process ICHORCNA_PLOT {
    tag      "${meta.id}"
    label    'process_low'

    // Copy-number figures drawn from ichorCNA, on the hg38 track, beside the
    // ichorkaryo JSON whose numbers they illustrate.
    //
    // This replaces MMPLOT on mmkaryo/mmbaf output. Until this change the
    // copy-number tab showed ichorkaryo's tables (ichorCNA, hg38, purity from
    // flow, ploidy pinned at 2) above figures drawn by mmkaryo (T2T, its own
    // fit) — two callers and two references on one page. One caller, one
    // reference, one tab.
    //
    // ichorcna_to_mmplot.py turns .cna.seg and .seg.txt into the positional
    // tables mmplot.py reads, carrying ichorCNA's integer copy number through
    // so mmplot draws the caller's states rather than re-deriving them from
    // log2. No allele data is supplied — ichorCNA is depth-only and allelic
    // state is measured at the eight wide panel windows by genebaf — so the
    // BAF panel is omitted rather than drawn empty.
    //
    // Sex is read from the ichorkaryo JSON (which took it from the sample
    // sheet) so the figure and the tables judge chrX/chrY against the same
    // constitutional expectation. The subtitle carries the ISCN string and the
    // ploidy class, the latter marked as computed rather than confirmed, per
    // the 15 September decision that ploidy is read from the plot.

    input:
    tuple val(meta), path(cna_seg), path(seg), path(karyotype_json), path(panel_bed)
    path cytobands
    path genes

    output:
    tuple val(meta), path("${meta.id}.cn_genome.png"),       emit: genome
    tuple val(meta), path("${meta.id}.cn_grid.png"),         emit: grid
    tuple val(meta), path("${meta.id}.cn_chr*.png"),         optional: true, emit: chromosomes
    tuple val(meta), path("${meta.id}.ichorcna.bins.tsv"),
                     path("${meta.id}.ichorcna.segments.tsv"), emit: tables
    path "versions.yml",                                     emit: versions

    script:
    def gene_arg   = genes.name != 'NO_FILE' ? "--genes ${genes}" : ''
    def target_arg = panel_bed.name != 'NO_FILE' ? "--targets ${panel_bed}" : ''
    def sex_fallback = meta.sex ?: 'XX'
    """
    set -euo pipefail

    ichorcna_to_mmplot.py \\
        --cna ${cna_seg} \\
        --seg ${seg} \\
        --out-prefix ${meta.id}.ichorcna

    # Sex and the subtitle come from the karyotype JSON so the figure cannot
    # disagree with the tables it sits above. The JSON's sex falls back to the
    # sample sheet, then to XX, matching ichorkaryo's own default.
    SEX=\$(python3 -c "import json; print(json.load(open('${karyotype_json}')).get('sex') or '${sex_fallback}')")
    SUBTITLE=\$(python3 - <<'PY'
import json
d = json.load(open("${karyotype_json}"))
parts = []
pc = d.get("ploidy_class")
if pc:
    tri = d.get("canonical_trisomies") or []
    tri = " ".join("+" + str(t).replace("chr", "") for t in tri)
    parts.append("%s%s (computed, to be confirmed)" % (pc, (" " + tri) if tri else ""))
iscn = d.get("ISCN_karyotype")
if iscn:
    parts.append(iscn)
print("; ".join(parts))
PY
)

    COMMON=(--no-baf --sex "\$SEX" --max-cn ${params.cn_plot_max_cn}
            --title "${meta.id}" --subtitle "\$SUBTITLE"
            --cytobands ${cytobands} ${gene_arg} ${target_arg})

    mmplot.py ${meta.id}.ichorcna.bins.tsv ${meta.id}.ichorcna.segments.tsv \\
        ${meta.id}.cn_genome.png "\${COMMON[@]}" --genome
    mmplot.py ${meta.id}.ichorcna.bins.tsv ${meta.id}.ichorcna.segments.tsv \\
        ${meta.id}.cn_grid.png "\${COMMON[@]}"

    # One page per chromosome. A chromosome with no segments makes mmplot exit
    # non-zero; that is skipped rather than rendered as an empty page.
    for C in ${params.cn_plot_chromosomes}; do
        mmplot.py ${meta.id}.ichorcna.bins.tsv ${meta.id}.ichorcna.segments.tsv \\
            ${meta.id}.cn_chr\${C}.png "\${COMMON[@]}" --chr chr\${C} || \\
            echo "  no data for chr\${C}, skipped" >&2
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | sed 's/Python //')
        matplotlib: \$(python3 -c "import matplotlib; print(matplotlib.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.cn_genome.png ${meta.id}.cn_grid.png \\
          ${meta.id}.ichorcna.bins.tsv ${meta.id}.ichorcna.segments.tsv
    echo '"${task.process}": stub' > versions.yml
    """
}
