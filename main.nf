#!/usr/bin/env nextflow

include { PREPARE_GENOMES } from './subworkflows/prepare_genomes'
include { ALIGN_AND_CHAIN  } from './subworkflows/align_and_chain'
include { LIFTOVER         } from './subworkflows/liftover'
include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow {
    if (params.help) {
        log.info """
        Usage:
          nextflow run ${workflow.projectDir} \\
            --id <chrom_pos.id> \\
            --ref_fa <original.fa> \\
            --query_fa <target.fa> \\
            --outdir <results_dir> \\
            -profile standard

        Required:
          --id        Input ID file, one chrom_pos ID per line
          --ref_fa    Original/reference FASTA
          --query_fa  Target/new FASTA

        Optional:
          --outdir            Output directory (default: results)
          --ref_fai           Existing index for --ref_fa
          --query_fai         Existing index for --query_fa
          --mapping           Two-column chromosome mapping TSV
          --pair_strategy     order or suffix (default: suffix)
          --split_threshold   Split chromosomes at this length (default: 100000000)
          --split_size        Sliding window size for large chromosomes (default: 10000000)
          --flank             snpcalling BED flank size (default: 100)
        """.stripIndent()
        System.exit(0)
    }

    validateParameters([parameters_schema: 'nextflow_schema.json'])
    log.info paramsSummaryLog([parameters_schema: 'nextflow_schema.json'], workflow)

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
