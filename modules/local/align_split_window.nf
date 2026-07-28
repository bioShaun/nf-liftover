nextflow.enable.types = true

record PairPaf {
    pair_id: String
    ref_chr: String
    query_chr: String
    paf: Path
}

process ALIGN_SPLIT_WINDOW {
    tag "${split_base}"
    label 'tool_aligner'
    label 'split_mem'
    stageInMode 'symlink'

    input:
    tuple(pair_id: String, ref_chr: String, query_chr: String, window_start: Integer, window_end: Integer, split_name: String, split_file: String, split_base: String, ref_chrom_fa: Path, ref_chrom_fai: Path, query_mmi: Path)

    output:
    paf: PairPaf = record(
        pair_id: pair_id,
        ref_chr: ref_chr,
        query_chr: query_chr,
        paf: file("${split_base}.paf")
    )
    versions: Path = file("versions.yml")

    script:
    def args = task.ext.args ?: ''
    def aligner = params.aligner ?: 'minimap2'
    """
    samtools faidx "${ref_chrom_fa}" "${ref_chr}:${window_start}-${window_end}" \\
      | awk -v name="${split_name}" 'NR == 1 { print ">" name; next } { print }' \\
      > "${split_file}"
    ${aligner} ${args} -t ${task.cpus} "${query_mmi}" "${split_file}" > "${split_base}.paf"

    {
      echo '"${task.process}":'
      echo "    samtools: \$(samtools --version | head -n 1 | awk '{ print \$2 }')"
      echo "    ${aligner}: \$(${aligner} --version)"
    } > versions.yml
    """

    stub:
    def aligner = params.aligner ?: 'minimap2'
    """
    touch "${split_base}.paf"
    {
      echo '"${task.process}":'
      echo "    samtools: 1.17"
      echo "    ${aligner}: stub"
    } > versions.yml
    """
}
