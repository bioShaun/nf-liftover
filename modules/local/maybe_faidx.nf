nextflow.enable.types = true

record IndexedGenome {
    role: String
    fasta: Path
    fai: Path
}

process MAYBE_FAIDX {
    tag "${role}"
    label 'tool_ngs'
    label 'small_mem'

    input:
    record(
        role: String,
        fasta: Path,
        fai: Path?
    )

    stage:
    // Fixed stage names: never interpolate user paths/basenames into the shell.
    stageAs fasta, 'input.fa'
    stageAs fai, 'provided.fai'

    output:
    indexed: IndexedGenome = record(role: role, fasta: file('genome.fa'), fai: file('genome.fa.fai'))
    versions: Path = file("versions.yml")

    script:
    if (fai) {
        """
        cp -L input.fa genome.fa
        cp -L provided.fai genome.fa.fai
        {
          echo '"${task.process}":'
          echo "    samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')"
          echo "    samtools_path: \$(command -v samtools)"
        } > versions.yml
        """
    } else {
        """
        cp -L input.fa genome.fa
        samtools faidx genome.fa
        {
          echo '"${task.process}":'
          echo "    samtools: \$(samtools --version | head -n 1 | awk '{print \$2}')"
          echo "    samtools_path: \$(command -v samtools)"
        } > versions.yml
        """
    }

    stub:
    """
    cp -L input.fa genome.fa 2>/dev/null || touch genome.fa
    if [ -f provided.fai ]; then
      cp -L provided.fai genome.fa.fai
    else
      awk '
        BEGIN { name=""; length_sum=0; seq_offset=0; line_bases=0; line_width=0 }
        /^>/ && name == "" {
          name=substr(\$0, 2)
          sub(/[ \\t].*/, "", name)
          seq_offset=length(\$0) + 1
          next
        }
        name != "" && \$0 !~ /^>/ {
          seq=\$0
          gsub(/[ \\t\\r]/, "", seq)
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
        }' genome.fa > genome.fa.fai
    fi
    {
      echo '"${task.process}":'
      echo "    samtools: 1.17"
      echo "    samtools_path: stub"
    } > versions.yml
    """
}
