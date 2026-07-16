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
}

workflow {
    missing_config = [
        pair_strategy: params.pair_strategy,
        align_mode: params.align_mode,
        split_threshold: params.split_threshold,
        split_size: params.split_size,
        flank: params.flank,
        conda_dir: params.conda_dir,
        conda_envs: params.conda_envs,
        max_cpus: params.max_cpus,
        max_memory: params.max_memory,
        max_time: params.max_time
    ].findAll { name, value -> value == null }.keySet()
    if (missing_config) {
        error "Missing pipeline configuration (${missing_config.join(', ')}). " +
            "Copy ${projectDir}/nextflow.config.example to ${projectDir}/nextflow.config " +
            "or provide an equivalent config with -c."
    }

    validateRequiredParams()

    input_paths = [
        id: params.id,
        vcf: params.vcf,
        ref_fa: params.ref_fa,
        query_fa: params.query_fa,
        ref_fai: params.ref_fai,
        query_fai: params.query_fai,
        mapping: params.mapping,
        split_bed: params.split_bed,
        split_genome_fai: params.split_genome_fai
    ].findAll { name, value -> value != null }
    missing_inputs = input_paths.findAll { name, value -> !file(value).exists() }.keySet()
    if (missing_inputs) {
        error "Input files do not exist: " + missing_inputs.collect { name ->
            "--${name} ${input_paths[name]}"
        }.join(', ')
    }

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
