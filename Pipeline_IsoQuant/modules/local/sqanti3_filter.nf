process SQANTI3_FILTER {
    tag "SQANTI3_Filter"
    label 'process_low'
    input:
    tuple val(group), path(original_classification)
    path corrected_fasta
    path corrected_gtf
    path filter_json

    output:
    path "sqanti_filter_out/${group}/*", emit: results
    path "versions.yml",                     emit: versions

    script:
    """

    mkdir -p sqanti_filter_out/${group}

    echo "1. Adding Roy score..."
    add_roy_score.py \\
        ${original_classification} \\
        ${group}_scored_classification.txt

    echo "2. Running SQANTI3 filter..."
    sqanti3_filter.py rules \\
        --sqanti_class ${group}_scored_classification.txt \\
        -j ${filter_json} \\
        --filter_isoforms ${corrected_fasta} \\
        --filter_gtf ${corrected_gtf} \\
        -d sqanti_filter_out/${group} \\
        -o ${group} \\
        --skip_report 

    SQ_VERSION=\$(sqanti3_filter.py --version 2>&1 | head -n1 | sed 's/SQANTI3 //g' | xargs || true)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$SQ_VERSION
    END_VERSIONS
    """
}
