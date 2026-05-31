process MINIMAP2 {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path genome
    val data_type // <--- ADD THIS

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix  = meta.id
    // Dynamic memory calculation for Samtools Sort
    def avail_mem = task.memory ? task.memory.toGiga() : 64
    def sort_mem  = Math.max(2, (avail_mem * 0.8 / task.cpus).intValue()) + 'G'
    
    def preset = data_type == 'pacbio_ccs' ? 'splice:hq' : 'splice -k14'

    """
    echo "DEBUG: Available Mem: ${avail_mem}GB, Threads: ${task.cpus}, Sort Buffer per thread: ${sort_mem}"

    # 1. Align
    # -t: threads
    # -a: output SAM
    # -x: dynamically set to splice:hq for PacBio or splice -k14 for Nanopore
    # -uf: force forward transcript strand (best for Iso-seq / cDNA)
    minimap2 \\
        -ax ${preset} \\
        -uf \\
        --secondary=no \\
        -MD \\
        -t ${task.cpus} \\
        ${genome} \\
        ${reads} \\
        > ${prefix}.sam

    # 2. Sort & Index
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${sort_mem} \\
        -o ${prefix}.sorted.bam \\
        ${prefix}.sam

    samtools index ${prefix}.sorted.bam

    # Cleanup intermediate SAM to save space
    rm ${prefix}.sam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}