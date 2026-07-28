nextflow.enable.types = true

process DERIVE_CHROM_PAIRS {
    tag 'chrom_pairs'
    label 'tool_py'
    label 'small_mem'

    input:
    ref_fai: Path
    query_fai: Path
    record(
        mapping_file: Path?
    )
    strategy: String

    stage:
    stageAs ref_fai, 'ref.fa.fai'
    stageAs query_fai, 'query.fa.fai'
    stageAs mapping_file, 'chrom_mapping.tsv'

    output:
    pairs: Path = file('chrom_pairs.tsv')
    versions: Path = file("versions.yml")

    script:
    def mapping_args = mapping_file ? "--mapping chrom_mapping.tsv" : ''
    """
    python ${projectDir}/bin/derive_chrom_pairs.py \\
      --ref-fai ref.fa.fai \\
      --query-fai query.fa.fai \\
      --output chrom_pairs.tsv \\
      --strategy "${strategy}" \\
      ${mapping_args}

    {
      echo '"${task.process}":'
      echo "    python: \$(python --version 2>&1 | awk '{ print \$2 }')"
    } > versions.yml
    """

    stub:
    """
    if [ -s chrom_mapping.tsv ]; then
      awk 'BEGIN { FS=OFS="\\t" } NF >= 2 && \$1 !~ /^#/ { print \$1, \$2; found=1; exit } END { if (!found) exit 1 }' chrom_mapping.tsv > chrom_pairs.tsv
    else
      awk 'BEGIN { FS=OFS="\\t" } NR == FNR { if (!ref) ref=\$1; next } { if (!query) query=\$1 } END { if (!ref || !query) exit 1; print ref, query }' ref.fa.fai query.fa.fai > chrom_pairs.tsv
    fi
    {
      echo '"${task.process}":'
      echo "    python: 3.12.0"
    } > versions.yml
    """
}
