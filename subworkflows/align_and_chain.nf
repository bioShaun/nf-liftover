nextflow.enable.types = true

process ALIGN_AND_CHAIN_PROCESS {
    tag 'align_and_chain'
    label 'tool_ngs'
    label 'large_mem'

    publishDir "${params.outdir}/chain", mode: 'copy', pattern: 'all.*'

    input:
    ref_fa: Path
    query_fa: Path
    ref_fai: Path
    query_fai: Path
    chrom_pairs: Path

    output:
    paf: Path = file('all.paf')
    chain: Path = file('all.chain')

    script:
    """
    mkdir -p chroms pafs

    while IFS=\$'\\t' read -r ref_chr query_chr; do
        [[ -z "\${ref_chr}" ]] && continue

        safe_ref=\$(printf '%s' "\${ref_chr}" | sed 's/[^A-Za-z0-9._-]/_/g')
        safe_query=\$(printf '%s' "\${query_chr}" | sed 's/[^A-Za-z0-9._-]/_/g')
        ref_chrom_fa="chroms/\${safe_ref}.ref.fa"
        query_chrom_fa="chroms/\${safe_query}.query.fa"

        samtools faidx "${ref_fa}" "\${ref_chr}" > "\${ref_chrom_fa}"
        samtools faidx "${query_fa}" "\${query_chr}" > "\${query_chrom_fa}"

        ref_len=\$(awk -v c="\${ref_chr}" '\$1 == c { print \$2; found=1; exit } END { if (!found) exit 1 }' "${ref_fai}")
        paf_prefix="pafs/\${safe_ref}_vs_\${safe_query}"

        if (( ref_len >= ${params.split_threshold} )); then
            seqkit sliding -W ${params.split_size} -s ${params.split_size} "\${ref_chrom_fa}" > "chroms/\${safe_ref}.split.fa"
            minimap2 -cx asm5 --cs -t ${task.cpus} "\${query_chrom_fa}" "chroms/\${safe_ref}.split.fa" > "\${paf_prefix}.split.paf"
            python ${projectDir}/bin/restore_split_paf.py "\${paf_prefix}.split.paf" "${ref_fai}" "\${paf_prefix}.paf"
        else
            minimap2 -cx asm5 --cs -t ${task.cpus} "\${query_chrom_fa}" "\${ref_chrom_fa}" > "\${paf_prefix}.paf"
        fi
    done < "${chrom_pairs}"

    cat pafs/*.paf > all.paf
    transanno minimap2chain --output all.chain all.paf
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
    alignment = ALIGN_AND_CHAIN_PROCESS(ref_fa, query_fa, ref_fai, query_fai, chrom_pairs)

    emit:
    paf   = alignment.paf
    chain = alignment.chain
}
