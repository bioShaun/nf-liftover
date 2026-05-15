nextflow.enable.types = true

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

    output:
    files = files('out/*')

    script:
    """
    mkdir -p out

    python ${projectDir}/bin/liftover_by_id.py \\
      "${id_file}" \\
      "${chain}" \\
      "${ref_fa}" \\
      "${query_fa}" \\
      out \\
      --flank ${params.flank}

    """
}

workflow LIFTOVER {
    take:
    id_file
    chain
    ref_fa
    query_fa
    query_fai

    main:
    lifted = LIFTOVER_BY_ID(id_file, chain, ref_fa, query_fa, query_fai)

    emit:
    files = lifted
}
