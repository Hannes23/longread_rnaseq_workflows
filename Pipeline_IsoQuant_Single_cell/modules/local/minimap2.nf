process MINIMAP2 {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path genome
    val data_type

    output:
    tuple val(meta), path("${meta.id}_chunk_${task.index}.sorted.bam"), path("${meta.id}_chunk_${task.index}.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    // Dynamic memory calculation for Samtools Sort
    def avail_mem = task.memory ? task.memory.toGiga() : 64
    def sort_mem  = Math.max(2, (avail_mem * 0.8 / task.cpus).intValue()) + 'G'
    
    // Dynamically assign the best minimap2 preset based on data type
    def preset = data_type == 'pacbio_ccs' ? 'splice:hq' : 'splice -k14'
    def prefix  = "${meta.id}_chunk_${task.index}"
    
    """
    set -euo pipefail  # ← Strict error handling
    
    echo "DEBUG: Available Mem: ${avail_mem}GB, Threads: ${task.cpus}, Sort Buffer per thread: ${sort_mem}"

    # 1. Align reads to genome
    minimap2 \\
        -ax ${preset} \\
        -uf \\
        --secondary=no \\
        -y \\
        -MD \\
        -t ${task.cpus} \\
        ${genome} \\
        ${reads} \\
        > ${prefix}.sam

    # 2. Convert SAM to sorted BAM immediately
    samtools sort \\
        -@ ${task.cpus} \\
        -m ${sort_mem} \\
        -o ${prefix}.sorted.bam \\
        ${prefix}.sam

    # 3. Create BAM index
    samtools index ${prefix}.sorted.bam

    # ← **CRITICAL**: Delete huge SAM file immediately (saves 10-50GB!)
    # This is the main optimization - SAM files are temporary and huge
    rm -f ${prefix}.sam

    # ← Verify outputs were created successfully
    if [ ! -f "${prefix}.sorted.bam" ]; then
        echo "ERROR: BAM file not created for ${prefix}!"
        exit 1
    fi
    
    if [ ! -f "${prefix}.sorted.bam.bai" ]; then
        echo "ERROR: BAM index not created for ${prefix}!"
        exit 1
    fi

    echo "Successfully created BAM for ${meta.id} chunk ${task.index}"

    # Generate version info
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version | head -n1 | sed 's/samtools //')
    END_VERSIONS
    """
}
