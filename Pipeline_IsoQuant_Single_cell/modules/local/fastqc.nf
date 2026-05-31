process FASTQC {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path "versions.yml"            , emit: versions

    script:
    """
     # A FASTQ record is 4 lines. 4,000,000 lines = 1,000,000 reads.
    zcat ${reads} | head -n 4000000 | gzip > subsampled_reads.fastq.gz

    fastqc \\
        -t ${task.cpus} \\
        subsampled_reads.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed -e "s/FastQC v//g" )
    END_VERSIONS
    """
}