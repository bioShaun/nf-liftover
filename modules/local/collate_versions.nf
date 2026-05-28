nextflow.enable.types = true

process COLLATE_VERSIONS {
    label 'tool_py'
    label 'small_mem'

    publishDir "${params.outdir}", mode: 'copy'

    input:
    nextflow_version: String
    versions_files

    output:
    versions_yml: Path = file('software_versions.yml')

    script:
    """
    cat <<-'EOF' > collate.py
import sys
lines = ['nextflow: "' + sys.argv[1] + '"']
seen_headers = set()
for fn in sys.argv[2:]:
    with open(fn) as f:
        content = f.read().strip()
        if not content:
            continue
        header = None
        for line in content.splitlines():
            if not line.strip():
                continue
            if line.startswith('"') or not line.startswith(' '):
                header = line.strip()
                if header in seen_headers:
                    header = None
                else:
                    seen_headers.add(header)
                    lines.append(line)
            elif header is not None:
                lines.append(line)
with open('software_versions.yml', 'w') as f:
    f.write('\\n'.join(lines) + '\\n')
EOF

    python3 collate.py "${nextflow_version}" ${versions_files.join(' ')}
    """

    stub:
    """
    touch software_versions.yml
    """
}
