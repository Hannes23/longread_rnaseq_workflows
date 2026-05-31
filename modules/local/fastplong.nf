process FASTPLONG {
    tag "$meta.id"                                          // create readabel label for log
    label 'process_medium'                                  

    input:
    tuple val(meta), path(reads)                            // standard nf-core input structure 

    output:
    tuple val(meta), path("*_fastplong.fastq"), emit: reads                     // bundles fastplong output with respective metadata, 'emit reads' --> creates named output to be used downstream
    tuple val(meta), path("*_failed_count.txt"), emit: count                    // do not output file of failed reads, only counts  
    path "versions.yml"                       , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"            // defines prefix to be used, if defined or not defined
    """
    gunzip -c $reads > input.fastq
    fastplong -i input.fastq -o ${prefix}_fastplong.fastq --failed_out ${prefix}_failed_reads.fastq -t $task.cpus
    # Check if file exists (it might not if 0 reads failed)
    if [ -f "${prefix}_failed_reads.fastq" ]; then
        # Count lines starting with @ (headers)
        count=\$(grep -c "^@" "${prefix}_failed_reads.fastq" || true)
        echo "Sample: $prefix" > "${prefix}_failed_count.txt"
        echo "Failed_Reads: \$count" >> "${prefix}_failed_count.txt"
        
        # DELETE the big file now that we have the number
        rm "${prefix}_failed_reads.fastq"
    else
       echo "Sample: $prefix" > "${prefix}_failed_count.txt"
    echo "Failed_Reads: 0" >> "${prefix}_failed_count.txt"
    fi

    rm input.fastq

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastplong: \$(fastplong --version 2>&1 | grep 'fastplong' | sed 's/fastplong //')
        grep: \$(grep --version | head -n 1)
    END_VERSIONS
    """
}