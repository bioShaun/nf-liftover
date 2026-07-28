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
from pathlib import Path

def clean_lines(text):
    # Normalize per-process versions snippets into valid YAML fragments.
    out = []
    for raw in text.splitlines():
        # Drop heredoc accidents (literal END) and blank lines
        stripped = raw.strip()
        if not stripped or stripped == "END":
            continue
        # Strip leading indentation from process headers produced by indented heredocs
        if stripped.startswith('"') and stripped.endswith(':'):
            out.append(stripped)
        elif raw.startswith(" ") or raw.startswith("\\t"):
            # Keep a consistent 4-space indent for nested keys
            key = stripped
            out.append("    " + key)
        else:
            out.append(stripped)
    return out

lines = ['nextflow: "' + sys.argv[1] + '"']
seen_headers = set()
for fn in sys.argv[2:]:
    path = Path(fn)
    if not path.is_file():
        continue
    content = path.read_text(encoding="utf-8", errors="replace")
    cleaned = clean_lines(content)
    if not cleaned:
        continue
    header = None
    for line in cleaned:
        if line.startswith('"') and line.endswith(':'):
            header = line
            if header in seen_headers:
                header = None
                continue
            seen_headers.add(header)
            lines.append(line)
        elif header is not None:
            lines.append(line)

Path("software_versions.yml").write_text("\\n".join(lines) + "\\n", encoding="utf-8")
EOF

    python3 collate.py "${nextflow_version}" ${versions_files.join(' ')}
    """

    stub:
    """
    echo 'nextflow: "stub"' > software_versions.yml
    """
}
