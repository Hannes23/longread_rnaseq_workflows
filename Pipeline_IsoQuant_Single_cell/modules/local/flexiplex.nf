process FLEXIPLEX_DISCOVER {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    val chemistry

    output:
    tuple val(meta), path("*_barcodes_counts.txt"), emit: counts
    path "versions.yml"                           , emit: versions

    script:
    """
    pigz -dc ${reads} -p ${task.cpus} | flexiplex -d ${chemistry} -f 0

    # Rename output to make it sample-specific
    mv flexiplex_barcodes_counts.txt ${meta.id}_barcodes_counts.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex -h | grep FLEXIPLEX | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}

process FLEXIPLEX_FILTER {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("*_barcodes_final.txt"), emit: final_barcodes
    path "versions.yml"                           , emit: versions

    script:
    """
    flexiplex-filter ${counts} > ${meta.id}_barcodes_final.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex -h | grep FLEXIPLEX | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}

process FLEXIPLEX_DEMUX {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads), path(barcodes)
    val chemistry
    val platform

    output:
    tuple val(meta), path("*.demultiplexed.fastq.gz"), emit: demux_fastq
    path "versions.yml", emit: versions

    script:
    // Platform-specific error tolerance
    def err_opts = platform.toUpperCase() == 'PACBIO' ? "-f 2 -e 1" : "-f 8 -e 2"

    // Robust decompression command
    def decompress_cmd = reads.name.endsWith('.gz') ? "pigz -dc ${reads}" : "cat ${reads}"
    def ext = params.compress ? "fastq.gz" : "fastq"
    def zip_cmd = params.compress ? "gzip" : "cat"

    """
    ${decompress_cmd} | flexiplex \
        -k ${barcodes} \
        ${err_opts} \
        -p ${task.cpus} | \
        ${zip_cmd} \
        > ${meta.id}.demultiplexed.${ext} \
        2> ${meta.id}_flexiplex_stats.txt

    # Only remove inputs if demuxed file was successfully created
    if [ -s "${meta.id}.demultiplexed.${ext}" ]; then
        rm -f "${reads}"
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flexiplex: \$(flexiplex -h | grep FLEXIPLEX | head -n 1 | awk '{print \$2}')
    END_VERSIONS
    """
}