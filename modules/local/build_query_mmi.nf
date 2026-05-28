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
    label 'tool_ngs'
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
    """
    minimap2 ${args} -t ${task.cpus} -d "${pair_fastas.pair_id}.query.mmi" "${pair_fastas.query_chrom_fa}"

    cat <<-END > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END
    """

    stub:
    """
    touch "${pair_fastas.pair_id}.query.mmi"
    cat <<-END > versions.yml
    "${task.process}":
        minimap2: 2.26
    END
    """
}
