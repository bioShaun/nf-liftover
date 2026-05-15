#!/usr/bin/env nextflow

include { PREPARE_GENOMES } from './subworkflows/prepare_genomes'
include { ALIGN_AND_CHAIN  } from './subworkflows/align_and_chain'
include { LIFTOVER         } from './subworkflows/liftover'

workflow {
    if (!params.id || !params.ref_fa || !params.query_fa) {
        error "Required parameters: --id, --ref_fa, --query_fa"
    }

    PREPARE_GENOMES(
        Channel.of(tuple('ref', params.ref_fa, file(params.ref_fa), params.ref_fai ?: '')),
        Channel.of(tuple('query', params.query_fa, file(params.query_fa), params.query_fai ?: '')),
        Channel.value(params.mapping ?: ''),
        Channel.value(params.pair_strategy)
    )

    ALIGN_AND_CHAIN(
        PREPARE_GENOMES.out.ref_fa,
        PREPARE_GENOMES.out.query_fa,
        PREPARE_GENOMES.out.ref_fai,
        PREPARE_GENOMES.out.query_fai,
        PREPARE_GENOMES.out.chrom_pairs
    )

    LIFTOVER(
        Channel.value(file(params.id)),
        ALIGN_AND_CHAIN.out.chain,
        PREPARE_GENOMES.out.ref_fa,
        PREPARE_GENOMES.out.query_fa,
        PREPARE_GENOMES.out.query_fai
    )
}

workflow.onComplete {
    log.info "nf-liftover complete: ${workflow.success ? 'OK' : 'FAILED'}"
}
