process SAMTOOLS_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("*.merged.sorted.bam"), path("*.merged.sorted.bam.bai"), emit: merged_bam
    path "versions.yml"                                                          , emit: versions

    script:
    def prefix = "${meta.id}"
    """
    # If there are multiple BAMs, merge them. If only one chunk exists, just copy it.
    if [ \$(echo ${bams} | wc -w) -gt 1 ]; then
        samtools merge -@ ${task.cpus} -o tmp.bam ${bams}
        samtools sort -@ ${task.cpus} -o ${prefix}.merged.sorted.bam tmp.bam
        rm tmp.bam
    else
        samtools sort -@ ${task.cpus} -o ${prefix}.merged.sorted.bam ${bams}
    fi

    # Create the .bai index file required by IsoQuant
    samtools index -@ ${task.cpus} ${prefix}.merged.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}