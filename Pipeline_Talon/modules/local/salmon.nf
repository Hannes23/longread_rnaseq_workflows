process SALMON_INDEX {
    label 'process_medium'

    input:
    path genome
    path transcripts

    output:
    path "salmon_index", emit: index
    path "versions.yml", emit: versions

    script:
    """
    grep "^>" $genome | cut -d " " -f 1 > decoys.txt
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcripts $genome > gentrome.fa
    salmon index -t gentrome.fa -d decoys.txt -p $task.cpus -i salmon_index --gencode
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}

process SALMON_QUANT {
    tag "$meta.id"
    label 'process_high'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    path "${meta.id}_quant"                   , emit: results
    path "${meta.id}_quant/logs/salmon_quant.log", emit: log
    path "versions.yml"                       , emit: versions

    script:
    """
    salmon quant -i $index -l A -r $reads --validateMappings -o ${meta.id}_quant -p $task.cpus
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        salmon: \$(salmon --version | sed 's/salmon //')
    END_VERSIONS
    """
}