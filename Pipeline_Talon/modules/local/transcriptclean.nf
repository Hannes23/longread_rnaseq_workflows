// Ensure CLEAN_GENOME is present
process CLEAN_GENOME {
    tag "$genome"
    label 'process_low' 
    input: path genome
    output: path "${genome.simpleName}.clean.fa", emit: clean_genome
    script:
    """
    awk '/^>/ {print \$1; next} {print}' $genome > ${genome.simpleName}.clean.fa
    """
}

process TRANSCRIPT_CLEAN {
    tag "${meta.id}.${chrom}"
    label 'process_high' // Splitting makes this safe for Medium resources
    
    input:
    tuple val(meta), val(chrom), path(bam) // Matches workflow tuple
    path clean_genome

    output:
    tuple val(meta), path("*_clean.sorted.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.id}_${chrom}"
    """
    # 1. Local Genome Copy
    cp -L $clean_genome local_genome.fa

    # 2. Run TranscriptClean (File Mode is fine for split files)
    samtools view -h $bam > input.sam

    TranscriptClean \
        --sam input.sam \
        --genome local_genome.fa \
        --outprefix ${prefix} \
        --threads $task.cpus \
        --canonOnly

    samtools sort -@ $task.cpus -o ${prefix}_clean.sorted.bam ${prefix}_clean.sam
    samtools index ${prefix}_clean.sorted.bam

    rm input.sam ${prefix}_clean.sam local_genome.fa*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: 2.0.2
    END_VERSIONS
    """
}