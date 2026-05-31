process ISOQUANT_DB {
    publishDir "${params.isoquant_db_dir}", mode: 'copy', saveAs: { filename -> filename }

    input:
    path gtf, stageAs: 'input.gtf'

    output:
    path "annotation.db", emit: db

    script:
    """
    if [ -s annotation.db ]; then
        echo "Using cached DB"
        exit 0
    fi

    python - <<EOF
import gffutils

gffutils.create_db(
    "input.gtf",
    "annotation.db",
    force=True,
    disable_infer_genes=True,
    disable_infer_transcripts=True,
    merge_strategy="create_unique"
)
EOF
    """
}
