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
    split_options: SplitLiftoverOptions

    output:
    files = files('out/*')

    script:
    def splitBedArg = split_options.split_bed ? "--split-bed \"${split_options.split_bed}\"" : ''
    def splitGenomeFaiArg = split_options.split_genome_fai ? "--split-genome-fai \"${split_options.split_genome_fai}\"" : ''
    """
    mkdir -p out

    python ${projectDir}/bin/liftover_by_id.py \\
      "${id_file}" \\
      "${chain}" \\
      "${ref_fa}" \\
      "${query_fa}" \\
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
