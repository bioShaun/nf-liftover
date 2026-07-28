nextflow.enable.types = true

process BUILD_ALIGN_PAIRS {
    tag 'align_pairs'
    label 'tool_ngs'
    label 'small_mem'

    input:
    chrom_pairs: Path
    ref_fai: Path

    stage:
    stageAs chrom_pairs, 'chrom_pairs.tsv'
    stageAs ref_fai, 'ref.fa.fai'

    output:
    pairs: Path = file('align_pairs.tsv')
    versions: Path = file("versions.yml")

    script:
    """
    awk 'BEGIN { FS=OFS="\\t" }
      NR == FNR { lengths[\$1]=\$2; next }
      NF == 0 || \$1 ~ /^#/ { next }
      {
        if (!(\$1 in lengths)) {
          printf "Chromosome %s not found in %s\\n", \$1, ARGV[1] > "/dev/stderr"
          exit 1
        }
        pair_id = \$1 "_vs_" \$2
        gsub(/[^A-Za-z0-9._-]/, "_", pair_id)
        print pair_id, \$1, \$2, lengths[\$1]
      }' ref.fa.fai chrom_pairs.tsv > align_pairs.tsv

    {
      echo '"${task.process}":'
      echo "    awk: \$(awk -W version 2>&1 | head -n 1 | awk '{ print \$3 }' || awk --version 2>&1 | head -n 1 | awk '{ print \$3 }')"
    } > versions.yml
    """

    stub:
    """
    awk 'BEGIN { FS=OFS="\\t" }
      NR == FNR { lengths[\$1]=\$2; next }
      NF == 0 || \$1 ~ /^#/ { next }
      {
        if (!(\$1 in lengths)) {
          printf "Chromosome %s not found in %s\\n", \$1, ARGV[1] > "/dev/stderr"
          exit 1
        }
        pair_id = \$1 "_vs_" \$2
        gsub(/[^A-Za-z0-9._-]/, "_", pair_id)
        print pair_id, \$1, \$2, lengths[\$1]
      }' ref.fa.fai chrom_pairs.tsv > align_pairs.tsv
    {
      echo '"${task.process}":'
      echo "    awk: 1.0.0"
    } > versions.yml
    """
}
