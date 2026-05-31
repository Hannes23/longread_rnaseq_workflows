process SAMTOOLS_MERGE {

    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*merged.bam"), path("*merged.bam.bai"), emit: bam


    script:
    """
    samtools merge -@ ${task.cpus} ${meta.id}.merged.bam ${bams.join(' ')}
    samtools index ${meta.id}.merged.bam
    """
}
