nextflow.enable.types = true

record IndexedGenome {
    role: String
    fasta: Path
    fai: Path
}

process MAYBE_FAIDX {
    tag "$role:${fasta.name}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    record(
        role: String,
        fasta: Path,
        fai: Path?
    )

    output:
    indexed: IndexedGenome = record(role: role, fasta: file(fasta.name), fai: file("${fasta.name}.fai"))
    versions: Path = file("versions.yml")

    script:
    if (fai) {
        """
        SRC_REAL=\$(readlink -f "${fai}")
        DST_REAL=\$(readlink -f "${fasta.name}.fai")
        if [ "\$SRC_REAL" != "\$DST_REAL" ]; then
            cp -L "${fai}" "${fasta.name}.fai"
        fi
        cat <<-END > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')
        END
        """
    } else {
        """
        samtools faidx "${fasta}"
        cat <<-END > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')
        END
        """
    }

    stub:
    """
    awk '
      BEGIN { name=""; length_sum=0; seq_offset=0; line_bases=0; line_width=0 }
      /^>/ && name == "" {
        name=substr(\$0, 2)
        sub(/[ \t].*/, "", name)
        seq_offset=length(\$0) + 1
        next
      }
      name != "" && \$0 !~ /^>/ {
        seq=\$0
        gsub(/[ \t\r]/, "", seq)
        if (line_bases == 0) {
          line_bases=length(seq)
          line_width=length(\$0) + 1
        }
        length_sum += length(seq)
      }
      END {
        if (name == "" || length_sum == 0) {
          exit 1
        }
        print name "\\t" length_sum "\\t" seq_offset "\\t" line_bases "\\t" line_width
      }' "${fasta}" > "${fasta.name}.fai"
    cat <<-END > versions.yml
    "${task.process}":
        samtools: 1.17
    END
    """
}
