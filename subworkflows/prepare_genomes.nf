nextflow.enable.types = true

record IndexedGenome {
    role: String
    fasta: Path
    fai: Path
}

process MAYBE_FAIDX {
    tag "$role:${fasta.name}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    record(
        role: String,
        fasta: Path,
        fai: Path?
    )

    output:
    indexed: IndexedGenome = record(role: role, fasta: file(fasta.name), fai: file("${fasta.name}.fai"))

    script:
    if (fai) {
        """
        cp "${fai}" "${fasta.name}.fai"
        """
    } else {
        """
        samtools faidx "${fasta}"
        """
    }
}

process DERIVE_CHROM_PAIRS {
    tag 'chrom_pairs'
    label 'tool_py'
    label 'small_mem'

    input:
    ref_fai: Path
    query_fai: Path
    record(
        mapping_file: Path?
    )
    strategy: String

    stage:
    stageAs ref_fai, 'ref.fa.fai'
    stageAs query_fai, 'query.fa.fai'
    stageAs mapping_file, 'chrom_mapping.tsv'

    output:
    pairs: Path = file('chrom_pairs.tsv')

    script:
    def mapping_args = mapping_file ? "--mapping chrom_mapping.tsv" : ''
    """
    python ${projectDir}/bin/derive_chrom_pairs.py \\
      --ref-fai ref.fa.fai \\
      --query-fai query.fa.fai \\
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
    indexed_genomes = MAYBE_FAIDX(ref_input.mix(query_input))

    ref_indexed = indexed_genomes.filter { genome -> genome.role == 'ref' }.first()
    query_indexed = indexed_genomes.filter { genome -> genome.role == 'query' }.first()

    pairs = DERIVE_CHROM_PAIRS(
        ref_indexed.map { genome -> genome.fai },
        query_indexed.map { genome -> genome.fai },
        mapping,
        strategy
    )

    emit:
    ref_fa      = ref_indexed.map { genome -> genome.fasta }
    ref_fai     = ref_indexed.map { genome -> genome.fai }
    query_fa    = query_indexed.map { genome -> genome.fasta }
    query_fai   = query_indexed.map { genome -> genome.fai }
    chrom_pairs = pairs
}
