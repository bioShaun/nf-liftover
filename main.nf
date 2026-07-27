#!/usr/bin/env nextflow

nextflow.enable.types = true

include { PREPARE_GENOMES } from './subworkflows/prepare_genomes'
include { ALIGN_AND_CHAIN  } from './subworkflows/align_and_chain'
include { LIFTOVER         } from './subworkflows/liftover'
include { LIFTOVER_VCF     } from './modules/local/liftover_vcf'
include { COLLATE_VERSIONS } from './modules/local/collate_versions'

def validateRequiredParams() {
    ['ref_fa', 'query_fa'].each { name ->
        if (!params[name]) {
            throw new IllegalArgumentException("Missing required parameter --${name}")
        }
    }
    if (!params.id && !params.vcf) {
        throw new IllegalArgumentException("Missing required input: provide --id or --vcf")
    }
    if (params.id && params.vcf) {
        throw new IllegalArgumentException("Provide only one input type: --id or --vcf")
    }
    def alignerEnvs = params.aligner_envs ?: [:]
    def aligner = params.aligner ?: 'minimap2'
    // Allow bare minimap2 even when aligner_envs is missing/incomplete (legacy configs).
    // Any other aligner must be explicitly mapped so conda beforeScript can locate its binary.
    if (aligner != 'minimap2' && !alignerEnvs.containsKey(aligner)) {
        throw new IllegalArgumentException(
            "Unsupported --aligner '${params.aligner}'; supported: minimap2, mm2plus.\n" +
            "If you upgraded an existing nextflow.config, copy the `aligner_envs` block from nextflow.config.example.")
    }
}

workflow {
    validateRequiredParams()

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

    if (params.vcf) {
        lifted = LIFTOVER_VCF(
            channel.value(file(params.vcf)),
            aligned.chain,
            prepared.ref_fa,
            prepared.query_fa,
            prepared.ref_fai,
            prepared.query_fai
        )
    } else {
        lifted = LIFTOVER(
            channel.value(file(params.id)),
            aligned.chain,
            prepared.ref_fa,
            prepared.query_fa,
            prepared.ref_fai,
            prepared.query_fai,
            split_liftover_ch
        )
    }

    ch_versions = prepared.versions
        .mix(aligned.versions)
        .mix(lifted.versions)
        .collect()

    COLLATE_VERSIONS(workflow.nextflow.version.toString(), ch_versions)
}
