nextflow.enable.types = true

process PAF_TO_CHAIN {
    tag 'all_chain'
    label 'tool_ngs'
    label 'small_mem'

    input:
    paf: Path

    output:
    chain: Path = file('all.chain')
    versions: Path = file("versions.yml")

    script:
    def args = task.ext.args ?: ''
    """
    transanno minimap2chain ${args} --output all.chain "${paf}"

    {
      echo '"${task.process}":'
      echo "    transanno: \$(transanno --version | awk '{ print \$2 }')"
    } > versions.yml
    """

    stub:
    """
    touch all.chain
    {
      echo '"${task.process}":'
      echo "    transanno: 1.2.0"
    } > versions.yml
    """
}
