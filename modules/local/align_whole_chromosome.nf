nextflow.enable.types = true

record PairFastas {
    pair_id: String
    ref_chr: String
    query_chr: String
    ref_len: Integer
    ref_chrom_fa: Path
    ref_chrom_fai: Path
    query_chrom_fa: Path
}

record PairPaf {
    pair_id: String
    ref_chr: String
    query_chr: String
    paf: Path
}

process ALIGN_WHOLE_CHROMOSOME {
    tag "${pair_fastas.pair_id}"
    label 'tool_aligner'
    label 'large_mem'

    input:
    pair_fastas: PairFastas

    output:
    paf: PairPaf = record(
        pair_id: pair_fastas.pair_id,
        ref_chr: pair_fastas.ref_chr,
        query_chr: pair_fastas.query_chr,
        paf: file("${pair_fastas.pair_id}.paf")
    )
    versions: Path = file("versions.yml")

    script:
    def args = task.ext.args ?: ''
    def aligner = params.aligner ?: 'minimap2'
    """
    ${aligner} ${args} -t ${task.cpus} "${pair_fastas.query_chrom_fa}" "${pair_fastas.ref_chrom_fa}" > "${pair_fastas.pair_id}.paf"

    {
      echo '"${task.process}":'
      echo "    ${aligner}: \$(${aligner} --version)"
      echo "    ${aligner}_path: \$(command -v ${aligner})"
    } > versions.yml
    """

    stub:
    def aligner = params.aligner ?: 'minimap2'
    """
    touch "${pair_fastas.pair_id}.paf"
    {
      echo '"${task.process}":'
      echo "    ${aligner}: stub"
      echo "    ${aligner}_path: stub"
    } > versions.yml
    """
}
