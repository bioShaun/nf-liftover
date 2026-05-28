#!/usr/bin/env nextflow

nextflow.enable.types = true

include { PREPARE_GENOMES } from './subworkflows/prepare_genomes'
include { ALIGN_AND_CHAIN  } from './subworkflows/align_and_chain'
include { LIFTOVER         } from './subworkflows/liftover'
include { COLLATE_VERSIONS } from './modules/local/collate_versions'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {
    validateParameters([parameters_schema: 'nextflow_schema.json'])
    log.info paramsSummaryLog([parameters_schema: 'nextflow_schema.json'], workflow)

    mapping_ch = params.mapping
        ? channel.value(record(mapping_file: file(params.mapping)))
        : channel.value(record(mapping_file: null))

    split_liftover_ch = channel.value(record(
        split_bed: params.split_bed ? file(params.split_bed) : null,
        split_genome_fai: params.split_genome_fai ? file(params.split_genome_fai) : null
    ))

    prepared = PREPARE_GENOMES(
        channel.of(record(role: 'ref', fasta: file(params.ref_fa), fai: params.ref_fai ? file(params.ref_fai) : null)),
        channel.of(record(role: 'query', fasta: file(params.query_fa), fai: params.query_fai ? file(params.query_fai) : null)),
        mapping_ch,
        channel.value(params.pair_strategy)
    )

    aligned = ALIGN_AND_CHAIN(
        prepared.ref_fa,
        prepared.query_fa,
        prepared.ref_fai,
        prepared.query_fai,
        prepared.chrom_pairs
    )

    lifted = LIFTOVER(
        channel.value(file(params.id)),
        aligned.chain,
        prepared.ref_fa,
        prepared.query_fa,
        prepared.query_fai,
        split_liftover_ch
    )

    ch_versions = prepared.versions
        .mix(aligned.versions)
        .mix(lifted.versions)
        .collect()

    COLLATE_VERSIONS(workflow.nextflow.version.toString(), ch_versions)
}
