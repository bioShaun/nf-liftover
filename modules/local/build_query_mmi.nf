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

record PairMmi {
    pair_id: String
    mmi: Path
}

process BUILD_QUERY_MMI {
    tag "${pair_fastas.pair_id}"
    label 'tool_aligner'
    label 'split_mem'

    input:
    pair_fastas: PairFastas

    output:
    pair_mmi: PairMmi = record(
        pair_id: pair_fastas.pair_id,
        mmi: file("${pair_fastas.pair_id}.query.mmi")
    )
    versions: Path = file("versions.yml")

    script:
    def args = task.ext.args ?: ''
    def aligner = params.aligner ?: 'minimap2'
    """
    ${aligner} ${args} -t ${task.cpus} -d "${pair_fastas.pair_id}.query.mmi" "${pair_fastas.query_chrom_fa}"

    cat <<-END > versions.yml
    "${task.process}":
        ${aligner}: \$(${aligner} --version)
    END
    """

    stub:
    def aligner = params.aligner ?: 'minimap2'
    """
    touch "${pair_fastas.pair_id}.query.mmi"
    cat <<-END > versions.yml
    "${task.process}":
        ${aligner}: stub
    END
    """
}
