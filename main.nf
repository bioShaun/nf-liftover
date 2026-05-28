#!/usr/bin/env nextflow

nextflow.enable.types = true

include { PREPARE_GENOMES } from './subworkflows/prepare_genomes'
include { ALIGN_AND_CHAIN  } from './subworkflows/align_and_chain'
include { LIFTOVER         } from './subworkflows/liftover'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

process COLLECT_SOFTWARE_VERSIONS {
    tag 'software_versions'
    label 'tool_py_ngs'
    label 'small_mem'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    nextflow_version: String

    output:
    versions: Path = file('software_versions.yml')

    script:
    """
    {
      echo 'nextflow: "${nextflow_version}"'
      echo "python: \\"\$(python --version 2>&1 | awk '{ print \$2 }')\\""
      echo "samtools: \\"\$(samtools --version | awk 'NR == 1 { print \$2 }')\\""
      echo "seqkit: \\"\$(seqkit version | awk '{ print \$2 }')\\""
      echo "minimap2: \\"\$(minimap2 --version)\\""
      echo "transanno: \\"\$(transanno --version | awk '{ print \$2 }')\\""
    } > software_versions.yml
    """
}

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

    LIFTOVER(
        channel.value(file(params.id)),
        aligned.chain,
        prepared.ref_fa,
        prepared.query_fa,
        prepared.query_fai,
        split_liftover_ch
    )

    COLLECT_SOFTWARE_VERSIONS(channel.value(workflow.nextflow.version.toString()))
}
