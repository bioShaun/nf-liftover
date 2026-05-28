nextflow.enable.types = true

include { BUILD_ALIGN_PAIRS } from '../modules/local/build_align_pairs'
include { EXTRACT_CHROM_FASTAS } from '../modules/local/extract_chrom_fastas'
include { BUILD_QUERY_MMI } from '../modules/local/build_query_mmi'
include { ALIGN_WHOLE_CHROMOSOME } from '../modules/local/align_whole_chromosome'
include { SPLIT_REF_CHROMOSOME } from '../modules/local/split_ref_chromosome'
include { ALIGN_SPLIT_WINDOW } from '../modules/local/align_split_window'
include { COMBINE_SPLIT_PAFS } from '../modules/local/combine_split_pafs'
include { COMBINE_ALL_PAFS } from '../modules/local/combine_all_pafs'
include { PAF_TO_CHAIN } from '../modules/local/paf_to_chain'

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

workflow ALIGN_AND_CHAIN {
    take:
    ref_fa
    query_fa
    ref_fai
    query_fai
    chrom_pairs

    main:
    align_pair_out = BUILD_ALIGN_PAIRS(chrom_pairs, ref_fai)
    align_pair_file = align_pair_out.pairs
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

    pair_fastas_out = EXTRACT_CHROM_FASTAS(ref_fa, query_fa, align_pairs)
    pair_fastas = pair_fastas_out.pair_fastas

    whole_pairs = pair_fastas.filter { pair ->
        params.align_mode == 'whole' || (params.align_mode == 'auto' && pair.ref_len < params.split_threshold)
    }
    split_pairs = pair_fastas.filter { pair ->
        params.align_mode == 'split' || (params.align_mode == 'auto' && pair.ref_len >= params.split_threshold)
    }

    whole_pafs_out = ALIGN_WHOLE_CHROMOSOME(whole_pairs)
    whole_pafs = whole_pafs_out.paf

    split_windows_out = SPLIT_REF_CHROMOSOME(split_pairs)
    split_windows_file = split_windows_out.plan
    split_windows = split_windows_file
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
    pair_query_mmis_out = BUILD_QUERY_MMI(split_pairs)
    pair_query_mmis = pair_query_mmis_out.pair_mmi
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
    split_align_outs = ALIGN_SPLIT_WINDOW(split_align_inputs)
    split_paf_groups = split_align_outs.paf
        .map { split_paf -> tuple(split_paf.pair_id, split_paf.ref_chr, split_paf.query_chr, split_paf.paf) }
        .groupTuple(by: [0, 1, 2])
    split_pafs_out = COMBINE_SPLIT_PAFS(split_paf_groups, ref_fai)
    split_pafs = split_pafs_out.paf

    collected_pafs = whole_pafs
        .mix(split_pafs)
        .map { pair_paf -> pair_paf.paf }
        .toSortedList { a, b -> a.name <=> b.name }
    all_paf_out = COMBINE_ALL_PAFS(collected_pafs)
    all_paf = all_paf_out.paf
    chain_out = PAF_TO_CHAIN(all_paf)
    chain = chain_out.chain

    versions = align_pair_out.versions
        .mix(pair_fastas_out.versions)
        .mix(whole_pafs_out.versions)
        .mix(split_windows_out.versions)
        .mix(pair_query_mmis_out.versions)
        .mix(split_align_outs.versions)
        .mix(split_pafs_out.versions)
        .mix(all_paf_out.versions)
        .mix(chain_out.versions)

    emit:
    paf      = all_paf
    chain    = chain
    versions = versions
}
