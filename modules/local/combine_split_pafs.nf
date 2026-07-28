nextflow.enable.types = true

record PairPaf {
    pair_id: String
    ref_chr: String
    query_chr: String
    paf: Path
}

process COMBINE_SPLIT_PAFS {
    tag "${pair_id}"
    label 'tool_py'
    label 'small_mem'

    input:
    tuple(pair_id: String, ref_chr: String, query_chr: String, split_pafs)
    ref_fai: Path

    output:
    paf: PairPaf = record(
        pair_id: pair_id,
        ref_chr: ref_chr,
        query_chr: query_chr,
        paf: file("${pair_id}.paf")
    )
    versions: Path = file("versions.yml")

    script:
    def splitPafList = split_pafs.collect { it.toString() }.sort().join('\n')
    """
    cat > split_pafs.list <<'EOF'
${splitPafList}
EOF
    xargs cat < split_pafs.list > "${pair_id}.split.paf"
    python ${projectDir}/bin/restore_split_paf.py "${pair_id}.split.paf" "${ref_fai}" "${pair_id}.paf"

    {
      echo '"${task.process}":'
      echo "    python: \$(python --version 2>&1 | awk '{ print \$2 }')"
    } > versions.yml
    """

    stub:
    """
    touch "${pair_id}.paf"
    {
      echo '"${task.process}":'
      echo "    python: 3.12.0"
    } > versions.yml
    """
}
