nextflow.enable.types = true

include { MAYBE_FAIDX } from '../modules/local/maybe_faidx'
include { DERIVE_CHROM_PAIRS } from '../modules/local/derive_chrom_pairs'

record IndexedGenome {
    role: String
    fasta: Path
    fai: Path
}

workflow PREPARE_GENOMES {
    take:
    ref_input
    query_input
    mapping
    strategy

    main:
    indexed_genomes_out = MAYBE_FAIDX(ref_input.mix(query_input))

    ref_indexed = indexed_genomes_out.indexed.filter { genome -> genome.role == 'ref' }.first()
    query_indexed = indexed_genomes_out.indexed.filter { genome -> genome.role == 'query' }.first()

    pairs_out = DERIVE_CHROM_PAIRS(
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
    chrom_pairs = pairs_out.pairs
    versions    = indexed_genomes_out.versions.mix(pairs_out.versions)
}
