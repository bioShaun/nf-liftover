nextflow.enable.types = true

record GenomeInput {
    role: String
    source_fasta: String
    fasta: Path
    fai_hint: String
}

record IndexedGenome {
    role: String
    fasta: Path
    fai: Path
}

record ChromMapping {
    has_mapping: Boolean
    mapping_file: Path?
}

process MAYBE_FAIDX {
    tag "$role:${fasta.name}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    record(
        role: String,
        source_fasta: String,
        fasta: Path,
        fai_hint: String
    )

    output:
    indexed: IndexedGenome = record(role: role, fasta: file(fasta.name), fai: file("${fasta.name}.fai"))

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
    ref_fai: Path
    query_fai: Path
    record(
        has_mapping: Boolean,
        mapping_file: Path?
    )
    strategy: String

    output:
    pairs: Path = file('chrom_pairs.tsv')

    script:
    def mapping_args = has_mapping ? "--mapping \"${mapping_file}\"" : ''
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

    ref_indexed = MAYBE_FAIDX.out.indexed.filter { genome -> genome.role == 'ref' }
    query_indexed = MAYBE_FAIDX.out.indexed.filter { genome -> genome.role == 'query' }

    DERIVE_CHROM_PAIRS(
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
    chrom_pairs = DERIVE_CHROM_PAIRS.out.pairs
}
