process SQANTI3_FILTER {
    tag "SQANTI3_Filter"
    

    input:
    tuple val(group_name), path(original_classification)
    path corrected_fasta
    path corrected_gtf
    path filter_json




    output:
    path "filter_output/*",                      emit: results
    path "versions.yml",                              emit: versions

    script:
    """
    set -euo pipefail

    mkdir -p filter_output

    echo "1. Adding Roy score..."
    add_roy_score.py \\
        ${original_classification} \\
        scored_classification.txt

    echo "2. Running SQANTI3 rules filter & Reporting..."
    sqanti3_filter.py rules \\
        --sqanti_class scored_classification.txt \\
        -j ${filter_json} \\
        --filter_isoforms ${corrected_fasta} \\
        --filter_gtf ${corrected_gtf} \\
        -d filter_output \\
        -o ${group_name} \\
        --skip_report 

    echo "3. Capturing version..."
    SQ_VERSION=\$(sqanti3_filter.py --version 2>&1 | head -n1 | sed 's/SQANTI3 //g' || echo "unknown")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$SQ_VERSION
    END_VERSIONS
    """
}