process LIFTOVER_BY_ID {
    tag "${id_file.baseName}"
    label 'tool_py_ngs'
    label 'medium_mem'

    publishDir "${params.outdir}/liftover", mode: 'copy'

    input:
    path id_file
    path chain
    path ref_fa
    path query_fa
    path query_fai

    output:
    path 'out/*', emit: files

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
    LIFTOVER_BY_ID(id_file, chain, ref_fa, query_fa, query_fai)

    emit:
    files = LIFTOVER_BY_ID.out.files
}
