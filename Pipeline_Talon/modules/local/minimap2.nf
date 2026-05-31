process MINIMAP2 {
    tag "$meta.id"

    input:
    tuple val(meta), path(reads)
    path genome

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix  = meta.id
    // Dynamic memory calculation for Samtools Sort
    // We allocate 80% of available memory to sorting, divided by the number of threads
    def avail_mem = task.memory ? task.memory.toGiga() : 64
    def sort_mem  = Math.max(2, (avail_mem * 0.8 / task.cpus).intValue()) + 'G'

    """
    echo "DEBUG: Available Mem: ${avail_mem}GB, Threads: ${task.cpus}, Sort Buffer per thread: ${sort_mem}"

    # 1. Align
    # -t: threads
    # -a: output SAM
    # -x splice: preset for splicing
    minimap2 \
        -ax splice \
        -uf \
        --secondary=no \
        -MD \
        -t ${task.cpus} \
        ${genome} \
        ${reads} \
        > ${prefix}.sam

    # 2. Sort & Index
    # -@: threads
    # -m: memory per thread (calculated dynamically above)
    samtools sort \
        -@ ${task.cpus} \
        -m ${sort_mem} \
        -o ${prefix}.sorted.bam \
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