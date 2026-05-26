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
    query_chrom_fa: Path
}

record SplitFastas {
    pair_id: String
    ref_chr: String
    query_chr: String
    query_chrom_fa: Path
    split_files: List<Path>
}

record SplitWindow {
    pair_id: String
    ref_chr: String
    query_chr: String
    query_chrom_fa: Path
    split_file: Path
}

record SplitPafGroup {
    pair_id: String
    ref_chr: String
    query_chr: String
    split_pafs: List<Path>
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
      }' "${ref_fai}" "${chrom_pairs}" > align_pairs.tsv
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

    output:
    pair_fastas: PairFastas = record(
        pair_id: pair.pair_id,
        ref_chr: pair.ref_chr,
        query_chr: pair.query_chr,
        ref_len: pair.ref_len,
        ref_chrom_fa: file("${pair.pair_id}.ref.fa"),
        query_chrom_fa: file("${pair.pair_id}.query.fa")
    )

    script:
    """
    samtools faidx "${ref_fa}" "${pair.ref_chr}" > "${pair.pair_id}.ref.fa"
    samtools faidx "${query_fa}" "${pair.query_chr}" > "${pair.pair_id}.query.fa"
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
    splits: SplitFastas = record(
        pair_id: pair_fastas.pair_id,
        ref_chr: pair_fastas.ref_chr,
        query_chr: pair_fastas.query_chr,
        query_chrom_fa: pair_fastas.query_chrom_fa,
        split_files: files('split_fa/*')
    )

    script:
    """
    mkdir -p split_fa
    seqkit sliding -W ${params.split_size} -s ${params.split_size} "${pair_fastas.ref_chrom_fa}" \\
      | seqkit split -i --by-id-prefix "" -O split_fa/
    """
}

process ALIGN_SPLIT_WINDOW {
    tag "${split_window.split_file.baseName}"
    label 'tool_ngs'
    label 'large_mem'

    input:
    split_window: SplitWindow

    output:
    paf: PairPaf = record(
        pair_id: split_window.pair_id,
        ref_chr: split_window.ref_chr,
        query_chr: split_window.query_chr,
        paf: file("${split_window.split_file.baseName}.paf")
    )

    script:
    """
    minimap2 -cx asm5 --cs -t ${task.cpus} "${split_window.query_chrom_fa}" "${split_window.split_file}" > "${split_window.split_file.baseName}.paf"
    """
}

process COMBINE_SPLIT_PAFS {
    tag "${split_group.pair_id}"
    label 'tool_py'
    label 'small_mem'

    input:
    split_group: SplitPafGroup
    ref_fai: Path

    output:
    paf: PairPaf = record(
        pair_id: split_group.pair_id,
        ref_chr: split_group.ref_chr,
        query_chr: split_group.query_chr,
        paf: file("${split_group.pair_id}.paf")
    )

    script:
    """
    cat ${split_group.split_pafs} > "${split_group.pair_id}.split.paf"
    python ${projectDir}/bin/restore_split_paf.py "${split_group.pair_id}.split.paf" "${ref_fai}" "${split_group.pair_id}.paf"
    """
}

process COMBINE_ALL_PAFS {
    tag 'all_paf'
    label 'tool_ngs'
    label 'small_mem'

    publishDir "${params.outdir}/chain", mode: 'copy', pattern: 'all.paf'

    input:
    paf_files: List<Path>

    output:
    paf: Path = file('all.paf')

    script:
    """
    cat ${paf_files} > all.paf
    """
}

process PAF_TO_CHAIN {
    tag 'all_chain'
    label 'tool_ngs'
    label 'large_mem'

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
        .flatMap { split ->
            split.split_files.collect { split_file ->
                record(
                    pair_id: split.pair_id,
                    ref_chr: split.ref_chr,
                    query_chr: split.query_chr,
                    query_chrom_fa: split.query_chrom_fa,
                    split_file: split_file
                )
            }
        }
    split_paf_groups = ALIGN_SPLIT_WINDOW(split_windows)
        .map { split_paf -> tuple(split_paf.pair_id, split_paf.ref_chr, split_paf.query_chr, split_paf.paf) }
        .groupTuple(by: [0, 1, 2])
        .map { pair_id, ref_chr, query_chr, split_pafs ->
            record(pair_id: pair_id, ref_chr: ref_chr, query_chr: query_chr, split_pafs: split_pafs)
        }
    split_pafs = COMBINE_SPLIT_PAFS(split_paf_groups, ref_fai)

    collected_pafs = whole_pafs
        .mix(split_pafs)
        .map { pair_paf -> pair_paf.paf }
        .collect()
    all_paf = COMBINE_ALL_PAFS(collected_pafs)
    chain = PAF_TO_CHAIN(all_paf)

    emit:
    paf   = all_paf
    chain = chain
}
