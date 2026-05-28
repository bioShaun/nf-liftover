nextflow.enable.types = true

record SplitLiftoverOptions {
    split_bed: Path?
    split_genome_fai: Path?
}

process LIFTOVER_BY_ID {
    tag "${id_file.baseName}"
    label 'tool_py_ngs'
    label 'medium_mem'

    publishDir "${params.outdir}/liftover", mode: 'copy', saveAs: { file -> file.tokenize('/').last() }

    input:
    id_file: Path
    chain: Path
    ref_fa: Path
    query_fa: Path
    query_fai: Path
    record(
        split_bed: Path?,
        split_genome_fai: Path?
    )

    stage:
    stageAs ref_fa, 'liftover_ref.fa'
    stageAs query_fa, 'liftover_query.fa'
    stageAs query_fai, 'liftover_query.fa.fai'
    stageAs split_bed, 'split_liftover.bed'
    stageAs split_genome_fai, 'split_liftover.genome.fai'

    output:
    files = files('out/*')

    script:
    def splitBedArg = split_bed ? "--split-bed split_liftover.bed" : ''
    def splitGenomeFaiArg = split_genome_fai ? "--split-genome-fai split_liftover.genome.fai" : ''
    """
    mkdir -p out

    python ${projectDir}/bin/liftover_by_id.py \\
      "${id_file}" \\
      "${chain}" \\
      liftover_ref.fa \\
      liftover_query.fa \\
      out \\
      --flank ${params.flank} \\
      ${splitBedArg} \\
      ${splitGenomeFaiArg}

    """
}

workflow LIFTOVER {
    take:
    id_file
    chain
    ref_fa
    query_fa
    query_fai
    split_options

    main:
    lifted = LIFTOVER_BY_ID(id_file, chain, ref_fa, query_fa, query_fai, split_options)

    emit:
    files = lifted
}
