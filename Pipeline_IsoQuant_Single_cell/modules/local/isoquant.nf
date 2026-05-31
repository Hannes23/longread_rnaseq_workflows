process ISOQUANT {
    tag "${group_name}" 
    label 'process_high_memory'

    input:
    tuple val(group_name), path(bams), path(bais), val(labels)
    path genome
    path gtf
    val data_type
    val umi_tag
    val barcode_tag
    val sc_mode

    output:
    tuple val(group_name), path("isoquant_out/${group_name}/${group_name}.extended_annotation.gtf"), emit: gtf
    tuple val(group_name), path("isoquant_out/${group_name}/*"), emit: full_output
    tuple val(group_name), path("isoquant_out/${group_name}/${group_name}.transcript_counts.tsv"),   emit: counts
    path "versions.yml", emit: versions

    script:
    def label_str = labels.join(' ')
    
    """
    python -m isoquant \\
        --reference ${genome} \\
        --genedb ${gtf} \\
        --bam ${bams} \\
        --labels ${label_str} \\
        --data_type ${data_type} \\
        --prefix ${group_name} \\
        --output isoquant_out \\
        --threads ${task.cpus} \\
        --complete_genedb \\
        --genedb_output . \\
        --sqanti_output \\
        --counts_format mtx \\
        --mode ${sc_mode} \\
        --barcoded_bam \\
        --barcode_tag ${barcode_tag} \\
        --umi_tag ${umi_tag} \\
        --transcript_quantification with_ambiguous

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        isoquant: \$(python -m isoquant --version | sed 's/IsoQuant //')
    END_VERSIONS
    """
}
