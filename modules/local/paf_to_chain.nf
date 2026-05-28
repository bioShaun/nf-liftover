nextflow.enable.types = true

process PAF_TO_CHAIN {
    tag 'all_chain'
    label 'tool_ngs'
    label 'small_mem'

    publishDir "${params.outdir}/chain", mode: 'copy', pattern: 'all.chain'

    input:
    paf: Path

    output:
    chain: Path = file('all.chain')
    versions: Path = file("versions.yml")

    script:
    def args = task.ext.args ?: ''
    """
    transanno minimap2chain ${args} --output all.chain "${paf}"

    cat <<-END > versions.yml
    "${task.process}":
        transanno: \$(transanno --version | awk '{ print \$2 }')
    END
    """

    stub:
    """
    touch all.chain
    cat <<-END > versions.yml
    "${task.process}":
        transanno: 1.2.0
    END
    """
}
