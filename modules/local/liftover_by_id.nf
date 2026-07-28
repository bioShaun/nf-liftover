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
    ref_fai: Path
    query_fai: Path
    record(
        split_bed: Path?,
        split_genome_fai: Path?
    )

    stage:
    // Keep the original id basename for output naming; only rename colliding FASTA/chain paths.
    stageAs chain, 'input.chain'
    stageAs ref_fa, 'liftover_ref.fa'
    stageAs query_fa, 'liftover_query.fa'
    stageAs ref_fai, 'liftover_ref.fa.fai'
    stageAs query_fai, 'liftover_query.fa.fai'
    stageAs split_bed, 'split_liftover.bed'
    stageAs split_genome_fai, 'split_liftover.genome.fai'

    output:
    files = files('out/*')
    versions: Path = file("versions.yml")

    script:
    def splitBedArg = split_bed ? "--split-bed split_liftover.bed" : ''
    def splitGenomeFaiArg = split_genome_fai ? "--split-genome-fai split_liftover.genome.fai" : ''
    def args = task.ext.args ?: ''
    def probeName = id_file.baseName
    """
    mkdir -p out

    python ${projectDir}/bin/liftover_by_id.py \\
      "${id_file}" \\
      input.chain \\
      liftover_ref.fa \\
      liftover_query.fa \\
      out \\
      --flank ${params.flank} \\
      --probe-name "${probeName}" \\
      ${splitBedArg} \\
      ${splitGenomeFaiArg} \\
      ${args}

    {
      echo '"${task.process}":'
      echo "    python: \$(python --version 2>&1 | awk '{ print \$2 }')"
      echo "    python_path: \$(command -v python)"
      echo "    transanno: \$(transanno --version | awk '{ print \$2 }')"
      echo "    transanno_path: \$(command -v transanno)"
    } > versions.yml
    """

    stub:
    """
    mkdir -p out
    touch "out/${id_file.baseName}.id"
    touch "out/${id_file.baseName}.bed"
    touch "out/${id_file.baseName}.pos.tsv"
    touch "out/${id_file.baseName}.snpcalling.bed"
    {
      echo '"${task.process}":'
      echo "    python: 3.12.0"
      echo "    transanno: 1.2.0"
    } > versions.yml
    """
}
