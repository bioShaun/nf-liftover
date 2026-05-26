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
          --mapping            Two-column chromosome mapping TSV
          --pair_strategy      order or suffix (default: suffix)
          --align_mode         auto, whole, or split (default: auto)
          --split_threshold    Split chromosomes at this length (default: 100000000)
          --split_size         Sliding window size for large chromosomes (default: 10000000)
          --split_bed          Optional split.bed for split-coordinate BED output
          --split_genome_fai   Optional split genome FASTA index for split BED sorting
          --flank              snpcalling BED flank size (default: 100)
        """.stripIndent()
        System.exit(0)
    }

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
        channel.of(record(role: 'ref', source_fasta: params.ref_fa, fasta: file(params.ref_fa), fai_hint: params.ref_fai ?: '')),
        channel.of(record(role: 'query', source_fasta: params.query_fa, fasta: file(params.query_fa), fai_hint: params.query_fai ?: '')),
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
