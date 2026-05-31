process ISOQUANT {

    tag "${group}"
    label 'process_high'

    input:
    tuple val(group), path(bam), path(bai), val(label)
    path genome
    path db 
    val data_type

    output:
    tuple val(group), path("isoquant_out/${group}/*.gtf"), emit: gtf
    tuple val(group), path("isoquant_out/${group}/*.tsv"), emit: counts
    path "versions.yml", emit: versions

    script:
    """
    python -m isoquant \
        --reference ${genome} \
        --genedb ${db} \
        --complete_genedb \
        --bam ${bam} \
        --labels ${label} \
        --data_type ${data_type} \
        --prefix ${group} \
        --output isoquant_out \
        --threads ${task.cpus} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(python -m isoquant --version | sed 's/IsoQuant //')
    END_VERSIONS
   """
}