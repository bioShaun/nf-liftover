nextflow.enable.types = true

record AlignPair {
    pair_id: String
    ref_chr: String
    query_chr: String
    ref_len: Integer
}

record PairFastas {
    pair_id: String
    ref_chr: String
    query_chr: String
    ref_len: Integer
    ref_chrom_fa: Path
    ref_chrom_fai: Path
    query_chrom_fa: Path
}

record PairMmi {
    pair_id: String
    mmi: Path
}

record PairPaf {
    pair_id: String
    ref_chr: String
    query_chr: String
    paf: Path
}

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
    """
}

process EXTRACT_CHROM_FASTAS {
    tag "${pair.pair_id}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    ref_fa: Path
    query_fa: Path
    pair: AlignPair

    stage:
    stageAs ref_fa, 'ref_genome.fa'
    stageAs query_fa, 'query_genome.fa'

    output:
    pair_fastas: PairFastas = record(
        pair_id: pair.pair_id,
        ref_chr: pair.ref_chr,
        query_chr: pair.query_chr,
        ref_len: pair.ref_len,
        ref_chrom_fa: file("${pair.pair_id}.ref.fa"),
        ref_chrom_fai: file("${pair.pair_id}.ref.fa.fai"),
        query_chrom_fa: file("${pair.pair_id}.query.fa")
    )

    script:
    """
    samtools faidx ref_genome.fa "${pair.ref_chr}" > "${pair.pair_id}.ref.fa"
    samtools faidx "${pair.pair_id}.ref.fa"
    samtools faidx query_genome.fa "${pair.query_chr}" > "${pair.pair_id}.query.fa"
    """
}

process BUILD_QUERY_MMI {
    tag "${pair_fastas.pair_id}"
    label 'tool_ngs'
    label 'split_mem'

    input:
    pair_fastas: PairFastas

    output:
    pair_mmi: PairMmi = record(
        pair_id: pair_fastas.pair_id,
        mmi: file("${pair_fastas.pair_id}.query.mmi")
    )

    script:
    """
    minimap2 -x asm5 -t ${task.cpus} -d "${pair_fastas.pair_id}.query.mmi" "${pair_fastas.query_chrom_fa}"
    """
}

process ALIGN_WHOLE_CHROMOSOME {
    tag "${pair_fastas.pair_id}"
    label 'tool_ngs'
    label 'large_mem'

    input:
    pair_fastas: PairFastas

    output:
    paf: PairPaf = record(
        pair_id: pair_fastas.pair_id,
        ref_chr: pair_fastas.ref_chr,
        query_chr: pair_fastas.query_chr,
        paf: file("${pair_fastas.pair_id}.paf")
    )

    script:
    """
    minimap2 -cx asm5 --cs -t ${task.cpus} "${pair_fastas.query_chrom_fa}" "${pair_fastas.ref_chrom_fa}" > "${pair_fastas.pair_id}.paf"
    """
}

process SPLIT_REF_CHROMOSOME {
    tag "${pair_fastas.pair_id}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    pair_fastas: PairFastas

    output:
    plan: Path = file('windows.tsv')

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
    """
}

process ALIGN_SPLIT_WINDOW {
    tag "${split_base}"
    label 'tool_ngs'
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

    script:
    """
    samtools faidx "${ref_chrom_fa}" "${ref_chr}:${window_start}-${window_end}" \\
      | awk -v name="${split_name}" 'NR == 1 { print ">" name; next } { print }' \\
      > "${split_file}"
    minimap2 -cx asm5 --cs -t ${task.cpus} "${query_mmi}" "${split_file}" > "${split_base}.paf"
    """
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

    script:
    def splitPafList = split_pafs.collect { it.toString() }.sort().join('\n')
    """
cat > split_pafs.list <<'EOF'
${splitPafList}
EOF
xargs cat < split_pafs.list > "${pair_id}.split.paf"
    python ${projectDir}/bin/restore_split_paf.py "${pair_id}.split.paf" "${ref_fai}" "${pair_id}.paf"
    """
}

process COMBINE_ALL_PAFS {
    tag 'all_paf'
    label 'tool_ngs'
    label 'small_mem'

    publishDir "${params.outdir}/chain", mode: 'copy', pattern: 'all.paf'

    input:
    paf_files

    output:
    paf: Path = file('all.paf')

    script:
    def pafList = paf_files.collect { it.toString() }.sort().join('\n')
    """
cat > pafs.list <<'EOF'
${pafList}
EOF
xargs cat < pafs.list > all.paf
    """
}

process PAF_TO_CHAIN {
    tag 'all_chain'
    label 'tool_ngs'
    label 'small_mem'

    publishDir "${params.outdir}/chain", mode: 'copy', pattern: 'all.chain'

    input:
    paf: Path

    output:
    chain: Path = file('all.chain')

    script:
    """
    transanno minimap2chain --output all.chain "${paf}"
    """
}

workflow ALIGN_AND_CHAIN {
    take:
    ref_fa
    query_fa
    ref_fai
    query_fai
    chrom_pairs

    main:
    align_pair_file = BUILD_ALIGN_PAIRS(chrom_pairs, ref_fai)
    align_pairs = align_pair_file
        .splitCsv(sep: '\t', header: false)
        .map { row ->
            record(
                pair_id: row[0] as String,
                ref_chr: row[1] as String,
                query_chr: row[2] as String,
                ref_len: row[3] as Integer
            )
        }

    pair_fastas = EXTRACT_CHROM_FASTAS(ref_fa, query_fa, align_pairs)

    whole_pairs = pair_fastas.filter { pair ->
        params.align_mode == 'whole' || (params.align_mode == 'auto' && pair.ref_len < params.split_threshold)
    }
    split_pairs = pair_fastas.filter { pair ->
        params.align_mode == 'split' || (params.align_mode == 'auto' && pair.ref_len >= params.split_threshold)
    }

    whole_pafs = ALIGN_WHOLE_CHROMOSOME(whole_pairs)

    split_windows = SPLIT_REF_CHROMOSOME(split_pairs)
        .splitCsv(sep: '\t', header: false)
        .map { row ->
            tuple(
                row[0] as String,
                row[1] as String,
                row[2] as String,
                row[3] as Integer,
                row[4] as Integer,
                row[5] as String,
                row[6] as String,
                row[7] as String
            )
        }
    pair_chrom_fastas = split_pairs.map { pf -> tuple(pf.pair_id, pf.ref_chrom_fa, pf.ref_chrom_fai) }
    pair_query_mmis = BUILD_QUERY_MMI(split_pairs)
        .map { pm -> tuple(pm.pair_id, pm.mmi) }
    split_align_inputs = split_windows
        .combine(pair_chrom_fastas, by: 0)
        .combine(pair_query_mmis, by: 0)
        .map { pair_id, ref_chr, query_chr, window_start, window_end, split_name, split_file, split_base, ref_chrom_fa, ref_chrom_fai, query_mmi ->
            tuple(
                pair_id,
                ref_chr,
                query_chr,
                window_start,
                window_end,
                split_name,
                split_file,
                split_base,
                ref_chrom_fa,
                ref_chrom_fai,
                query_mmi
            )
        }
    split_paf_groups = ALIGN_SPLIT_WINDOW(split_align_inputs)
        .map { split_paf -> tuple(split_paf.pair_id, split_paf.ref_chr, split_paf.query_chr, split_paf.paf) }
        .groupTuple(by: [0, 1, 2])
    split_pafs = COMBINE_SPLIT_PAFS(split_paf_groups, ref_fai)

    collected_pafs = whole_pafs
        .mix(split_pafs)
        .map { pair_paf -> pair_paf.paf }
        .toSortedList { a, b -> a.name <=> b.name }
    all_paf = COMBINE_ALL_PAFS(collected_pafs)
    chain = PAF_TO_CHAIN(all_paf)

    emit:
    paf   = all_paf
    chain = chain
}
