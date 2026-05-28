nextflow.enable.types = true

process COMBINE_ALL_PAFS {
    tag 'all_paf'
    label 'tool_ngs'
    label 'small_mem'

    publishDir "${params.outdir}/chain", mode: 'copy', pattern: 'all.paf'

    input:
    paf_files

    output:
    paf: Path = file('all.paf')
    versions: Path = file("versions.yml")

    script:
    def pafList = paf_files.collect { it.toString() }.sort().join('\n')
    """
    cat > pafs.list <<'EOF'
${pafList}
EOF
    xargs cat < pafs.list > all.paf

    cat <<-END > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | awk '{ print \$4 }' || echo "5.0")
    END
    """

    stub:
    """
    touch all.paf
    cat <<-END > versions.yml
    "${task.process}":
        bash: 5.2
    END
    """
}
