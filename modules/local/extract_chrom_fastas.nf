nextflow.enable.types = true

record AlignPair {
    pair_id: String
    ref_chr: String
    query_chr: String
    ref_len: Integer
}

record PairFastas {
    pair_id: String
    ref_chr: String
    query_chr: String
    ref_len: Integer
    ref_chrom_fa: Path
    ref_chrom_fai: Path
    query_chrom_fa: Path
}

process EXTRACT_CHROM_FASTAS {
    tag "${pair.pair_id}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    ref_fa: Path
    query_fa: Path
    ref_fai: Path
    query_fai: Path
    pair: AlignPair

    stage:
    stageAs ref_fa, 'ref_genome.fa'
    stageAs query_fa, 'query_genome.fa'
    stageAs ref_fai, 'ref_genome.fa.fai'
    stageAs query_fai, 'query_genome.fa.fai'

    output:
    pair_fastas: PairFastas = record(
        pair_id: pair.pair_id,
        ref_chr: pair.ref_chr,
        query_chr: pair.query_chr,
        ref_len: pair.ref_len,
        ref_chrom_fa: file("${pair.pair_id}.ref.fa"),
        ref_chrom_fai: file("${pair.pair_id}.ref.fa.fai"),
        query_chrom_fa: file("${pair.pair_id}.query.fa")
    )
    versions: Path = file("versions.yml")

    script:
    """
    samtools faidx ref_genome.fa "${pair.ref_chr}" > "${pair.pair_id}.ref.fa"
    samtools faidx "${pair.pair_id}.ref.fa"
    samtools faidx query_genome.fa "${pair.query_chr}" > "${pair.pair_id}.query.fa"

    {
      echo '"${task.process}":'
      echo "    samtools: \$(samtools --version | head -n 1 | awk '{ print \$2 }')"
    } > versions.yml
    """

    stub:
    """
    touch "${pair.pair_id}.ref.fa"
    touch "${pair.pair_id}.ref.fa.fai"
    touch "${pair.pair_id}.query.fa"
    {
      echo '"${task.process}":'
      echo "    samtools: 1.17"
    } > versions.yml
    """
}
