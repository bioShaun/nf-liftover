process MAYBE_FAIDX {
    tag "$role:${fasta.name}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    tuple val(role), val(source_fasta), path(fasta), val(fai_hint)

    output:
    tuple val(role), path(fasta), path("${fasta.name}.fai"), emit: indexed

    script:
    """
    source_fai="${fai_hint}"
    if [[ -z "\${source_fai}" ]]; then
        source_fai="${source_fasta}.fai"
    fi

    if [[ -s "\${source_fai}" ]]; then
        cp "\${source_fai}" "${fasta.name}.fai"
    else
        samtools faidx "${fasta}"
    fi
    """
}

process DERIVE_CHROM_PAIRS {
    tag 'chrom_pairs'
    label 'tool_py'
    label 'small_mem'

    input:
    path ref_fai
    path query_fai
    val mapping
    val strategy

    output:
    path 'chrom_pairs.tsv', emit: pairs

    script:
    def mapping_args = mapping ? "--mapping ${mapping}" : ''
    """
    python ${projectDir}/bin/derive_chrom_pairs.py \\
      --ref-fai "${ref_fai}" \\
      --query-fai "${query_fai}" \\
      --output chrom_pairs.tsv \\
      --strategy "${strategy}" \\
      ${mapping_args}
    """
}

workflow PREPARE_GENOMES {
    take:
    ref_input
    query_input
    mapping
    strategy

    main:
    MAYBE_FAIDX(ref_input.mix(query_input))

    ref_indexed = MAYBE_FAIDX.out.indexed.filter { role, fasta, fai -> role == 'ref' }
    query_indexed = MAYBE_FAIDX.out.indexed.filter { role, fasta, fai -> role == 'query' }

    DERIVE_CHROM_PAIRS(
        ref_indexed.map { role, fasta, fai -> fai },
        query_indexed.map { role, fasta, fai -> fai },
        mapping,
        strategy
    )

    emit:
    ref_fa      = ref_indexed.map { role, fasta, fai -> fasta }
    ref_fai     = ref_indexed.map { role, fasta, fai -> fai }
    query_fa    = query_indexed.map { role, fasta, fai -> fasta }
    query_fai   = query_indexed.map { role, fasta, fai -> fai }
    chrom_pairs = DERIVE_CHROM_PAIRS.out.pairs
}
