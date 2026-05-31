process SQANTI3 {
    tag "SQANTI3_QC"
    label 'process_high'
    
    input:
    tuple val(group), path(transcript_gtf), path(abundance)
    path ref_gtf
    path genome

    output:
    path("sqanti_out/${group}/*"),                    emit: results
    tuple val(group), path("sqanti_out/${group}/*_classification.txt"), emit: original_classification
    path("sqanti_out/${group}/*_corrected.fasta"),    emit: fasta
    path("sqanti_out/${group}/*_corrected.gtf"),      emit: corrected_gtf
    path("sqanti_out/${group}/*_junctions.txt"),      emit: junctions
    path "versions.yml",                                                        emit: versions

    script:
    """
    mkdir -p sqanti_out/${group}

    echo "1. Filtering GTF..."
    grep -E "^chr[0-9XYM]+[ \\t]" ${transcript_gtf} > filtered_isoquant.gtf || true
    grep -E "^chr[0-9XYM]+[ \\t]" ${ref_gtf} > filtered_ref.gtf || true

    if [ ! -s filtered_isoquant.gtf ]; then cp ${transcript_gtf} filtered_isoquant.gtf; fi
    if [ ! -s filtered_ref.gtf ]; then cp ${ref_gtf} filtered_ref.gtf; fi

    echo "2. Reformatting IsoQuant Abundance..."
    awk -F'\\t' '
    BEGIN { OFS="," }
    NR==1 {
        printf "id"
        for(i=12; i<=NF; i++) { printf "," \$i }
        printf "\\n"
    }
    NR>1 {
        id=\$4
        if(id=="None" || id=="" || id=="nan") id=\$1
        printf id
        for(i=12; i<=NF; i++) { printf "," \$i }
        printf "\\n"
    }' ${abundance} > sqanti_abundance.csv

    echo "3. Running SQANTI3 QC..."
    sqanti3_qc.py \\
        --isoforms filtered_isoquant.gtf \\
        --refGTF filtered_ref.gtf \\
        --refFasta ${genome} \\
        -fl sqanti_abundance.csv \\
        -d sqanti_out/${group} \\
        -o ${group}_sqanti \\
        -t ${task.cpus} \\
        --report skip

    SQ_VERSION=\$(sqanti3_qc.py --version 2>&1 | head -n1 | sed 's/SQANTI3 //')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \${SQ_VERSION}
    END_VERSIONS
    """
}
