process SPLIT_BAM_BY_CHR {
    tag "$meta.id"
    label 'process_medium'

   input:
    tuple val(meta), path(bam), path(bai) // Ensure input matches MINIMAP2 output

    output:
    tuple val(meta), path("${meta.id}_*.bam"), emit: chromosome_bams

    script:
    """
    # Use existing index
    # Get list of chromosomes
    samtools idxstats $bam | cut -f 1 | grep -E '^(chr)?[0-9XYM]+\$' > chrom_list.txt

    while read chrom; do
        [ -z "\$chrom" ] && continue
        
        # Check if chromosome has reads
        count=\$(samtools view -c $bam "\$chrom")
        
        if [ "\$count" -gt 0 ]; then
            echo "Splitting \$chrom..."
            samtools view -b $bam "\$chrom" > "${meta.id}_\${chrom}.bam"
        fi
    done < chrom_list.txt
    """
}