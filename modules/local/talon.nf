process TALON_INIT {
    input:
    path gtf

    output:
    path "talon.db", emit: db
    path "versions.yml", emit: versions

    script:
    """
    echo "Forcing a fresh DB rebuild to clear poisoned cache - attempt 2!"

    talon_initialize_database \\
        --f $gtf \\
        --g hg38 \\
        --a gencode_v49 \\
        --o talon

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        talon: \$(talon --version 2>&1 | sed 's/TALON //')
    END_VERSIONS
    """
}

process TALON_RUN {
    input:
    path bam
    path db
    path genome
    val platform

    output:
    path "talon_run.db", emit: db   
    path "versions.yml", emit: versions

    script:
    """
    # 1. Make a PHYSICAL copy of the database so we don't corrupt the symlink
    cp $db talon_run.db              

    # 2. Generate config
    # We escape \$ for Bash variables, but keep \${} for Nextflow variables
    BAM_FILE="${bam}"
    s_id="\$(basename "\$BAM_FILE" .merged.bam)_\$(date +%F)"
    
    echo "\$s_id,\$s_id,SequelII,\$BAM_FILE" > talon_config.csv

    # 3. Run TALON on the copied database
    talon \\
        --f talon_config.csv \\
        --db talon_run.db \\
        --build hg38 \\
        --o talon_results \\
        --threads ${task.cpus}       
     
    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        talon: \$(talon --version 2>&1 | sed 's/TALON //')
    END_VERSIONS
    """
}

process TALON_CREATE_GTF {
    input:
    path db 

    output:
    path "talon_filtered_talon.gtf", emit: gtf
    path "talon_filtered_talon_abundance_filtered.tsv", emit: abundance
    path "versions.yml", emit: versions

    script:
    """
    echo "1. Filtering Transcripts..."
    talon_filter_transcripts \\
        --db $db \\
        --annot gencode_v49 \\
        --minCount 5 \\
        --o talon_whitelist.csv

    echo "2. Creating Filtered GTF..."
    talon_create_GTF \\
        --db $db \\
        --build hg38 \\
        --annot gencode_v49 \\
        --whitelist talon_whitelist.csv \\
        --o talon_filtered

    echo "3. Creating Filtered Abundance..."
    talon_abundance \\
        --db $db \\
        --annot gencode_v49 \\
        --build hg38 \\
        --whitelist talon_whitelist.csv \\
        --o talon_filtered

    cat <<- END_VERSIONS > versions.yml
    "${task.process}":
        talon: \$(talon --version 2>&1 | sed 's/TALON //')
    END_VERSIONS
    """
}