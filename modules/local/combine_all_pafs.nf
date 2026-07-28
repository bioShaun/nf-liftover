nextflow.enable.types = true

process COMBINE_ALL_PAFS {
    tag 'all_paf'
    label 'tool_ngs'
    label 'small_mem'

    // Publish merged PAF only when --publish_paf is true (default false for production).
    publishDir "${params.outdir}/chain", mode: 'copy', pattern: 'all.paf', enabled: params.publish_paf

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

    {
      echo '"${task.process}":'
      echo "    bash: \$(bash --version | head -n 1 | awk '{ print \$4 }' || echo 5.0)"
    } > versions.yml
    """

    stub:
    """
    touch all.paf
    {
      echo '"${task.process}":'
      echo "    bash: 5.2"
    } > versions.yml
    """
}
