nextflow.enable.types = true

record PairFastas {
    pair_id: String
    ref_chr: String
    query_chr: String
    ref_len: Integer
    ref_chrom_fa: Path
    ref_chrom_fai: Path
    query_chrom_fa: Path
}

process SPLIT_REF_CHROMOSOME {
    tag "${pair_fastas.pair_id}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    pair_fastas: PairFastas

    output:
    plan: Path = file('windows.tsv')
    versions: Path = file("versions.yml")

    script:
    """
    : > windows.tsv
    start=1
    while (( start <= ${pair_fastas.ref_len} )); do
        end=\$(( start + ${params.split_size} - 1 ))
        if (( end > ${pair_fastas.ref_len} )); then
            end=${pair_fastas.ref_len}
        fi
        split_name="${pair_fastas.ref_chr}_sliding:\${start}-\${end}"
        split_base=\$(printf '%s' "\${split_name}" | sed 's/[^A-Za-z0-9._-]/_/g')
        printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
          "${pair_fastas.pair_id}" \\
          "${pair_fastas.ref_chr}" \\
          "${pair_fastas.query_chr}" \\
          "\${start}" \\
          "\${end}" \\
          "\${split_name}" \\
          "\${split_base}.fa" \\
          "\${split_base}" >> windows.tsv
        start=\$(( start + ${params.split_size} ))
    done

    cat <<-END > versions.yml
    "${task.process}":
        bash: \$(bash --version | head -n 1 | awk '{ print \$4 }' || echo "5.0")
    END
    """

    stub:
    """
    end=${params.split_size}
    if (( end > ${pair_fastas.ref_len} )); then
        end=${pair_fastas.ref_len}
    fi
    split_name="${pair_fastas.ref_chr}_sliding:1-\${end}"
    split_base=\$(printf '%s' "\${split_name}" | sed 's/[^A-Za-z0-9._-]/_/g')
    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \\
      "${pair_fastas.pair_id}" \\
      "${pair_fastas.ref_chr}" \\
      "${pair_fastas.query_chr}" \\
      "1" \\
      "\${end}" \\
      "\${split_name}" \\
      "\${split_base}.fa" \\
      "\${split_base}" > windows.tsv
    cat <<-END > versions.yml
    "${task.process}":
        bash: 5.2
    END
    """
}
