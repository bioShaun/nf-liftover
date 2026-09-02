nextflow.enable.types = true

process LIFTOVER_VCF {
    tag "${vcf_file.baseName}"
    label 'tool_ngs'
    label 'medium_mem'

    input:
    vcf_file: Path
    chain: Path
    ref_fa: Path
    query_fa: Path
    ref_fai: Path
    query_fai: Path

    stage:
    stageAs vcf_file, 'input.vcf'
    stageAs chain, 'input.chain'
    stageAs ref_fa, 'liftover_ref.fa'
    stageAs query_fa, 'liftover_query.fa'
    stageAs ref_fai, 'liftover_ref.fa.fai'
    stageAs query_fai, 'liftover_query.fa.fai'

    output:
    files = files('out/*')
    versions: Path = file("versions.yml")

    script:
    def args = task.ext.args ?: ''
    """
    mkdir -p out

    prefix=\$(basename "${vcf_file}")
    prefix=\${prefix%.vcf.gz}
    prefix=\${prefix%.vcf}
    prefix=\${prefix%.bcf}
    out_vcf="out/\${prefix}.vcf.gz"
    rejected_vcf="out/rejected.\${prefix}.vcf.gz"
    sorted_out_vcf="out/\${prefix}.sorted.tmp.vcf.gz"
    sorted_rejected_vcf="out/rejected.\${prefix}.sorted.tmp.vcf.gz"

    transanno liftvcf \
      --original-assembly liftover_ref.fa \
      --new-assembly liftover_query.fa \
      --chain "${chain}" \
      --vcf "${vcf_file}" \
      --output "\${out_vcf}" \
      --fail "\${rejected_vcf}" \
      ${args}

    if [ -s "\${out_vcf}" ]; then
      bcftools sort -Oz -o "\${sorted_out_vcf}" -T "\${TMPDIR}/bcftools-sort-success.XXXXXX" "\${out_vcf}"
      mv -f "\${sorted_out_vcf}" "\${out_vcf}"
      tabix -f -p vcf "\${out_vcf}"
    fi
    if [ -s "\${rejected_vcf}" ]; then
      bcftools sort -Oz -o "\${sorted_rejected_vcf}" -T "\${TMPDIR}/bcftools-sort-rejected.XXXXXX" "\${rejected_vcf}"
      mv -f "\${sorted_rejected_vcf}" "\${rejected_vcf}"
      tabix -f -p vcf "\${rejected_vcf}"
    fi

    {
      echo '"${task.process}":'
      echo "    transanno: \$(transanno --version | awk '{ print \$2 }')"
      echo "    tabix: \$(tabix --version 2>&1 | head -n 1 | awk '{ print \$2 }')"
    } > versions.yml
    """

    stub:
    """
    mkdir -p out
    prefix=\$(basename "${vcf_file}")
    prefix=\${prefix%.vcf.gz}
    prefix=\${prefix%.vcf}
    prefix=\${prefix%.bcf}
    touch "out/\${prefix}.vcf.gz"
    touch "out/\${prefix}.vcf.gz.tbi"
    touch "out/rejected.\${prefix}.vcf.gz"
    touch "out/rejected.\${prefix}.vcf.gz.tbi"
    {
      echo '"${task.process}":'
      echo "    transanno: 1.2.0"
      echo "    tabix: 1.17"
    } > versions.yml
    """
}
