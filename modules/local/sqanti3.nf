process SQANTI3 {
    input:
    path transcript_gtf
    path abundance
    path ref_gtf
    path genome

    output:
    path "sqanti_out/sqanti_output*",            emit: results
    path "sqanti_out/*_classification.txt",      emit: original_classification
    path "sqanti_out/*_corrected.fasta",         emit: fasta
    path "sqanti_out/*_corrected.gtf",           emit: corrected_gtf
    path "sqanti_out/*_junctions.txt",           emit: junctions   
    path "sqanti_out/*.params.txt",              emit: sqanti_params // FIXED: Renamed to avoid Nextflow keyword collision
    path "versions.yml",                         emit: versions  

    script:
    def sqanti_dir = "/data/gent/vo/000/gvo00082/vsc48905/Nextflow/apps/SQANTI3"

    """
    # 0. Set up safe environment to bypass broken container defaults
    export PYTHONPATH="${sqanti_dir}/cDNA_Cupcake/sequence/:\${PYTHONPATH:-}"
    export PATH="/env/SQANTI3/bin:\${PATH:-}"
    export LD_LIBRARY_PATH="/env/SQANTI3/lib:\${LD_LIBRARY_PATH:-}"

    mkdir -p sqanti_out

    echo "1. Filtering GTF..."
    grep -E "^chr[0-9XYM]+[ \\t]" ${transcript_gtf} > filtered_talon.gtf || true
    grep -E "^chr[0-9XYM]+[ \\t]" ${ref_gtf} > filtered_ref.gtf || true

    if [ ! -s filtered_talon.gtf ]; then cp ${transcript_gtf} filtered_talon.gtf; fi
    if [ ! -s filtered_ref.gtf ]; then cp ${ref_gtf} filtered_ref.gtf; fi

    echo "2. Reformatting Talon Abundance..."
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
    }' $abundance > sqanti_abundance.csv

    echo "3. Running SQANTI3 QC..."
    # Force the use of the host python and script
    /env/SQANTI3/bin/python ${sqanti_dir}/sqanti3_qc.py \\
        filtered_talon.gtf \\
        filtered_ref.gtf \\
        $genome \\
        --fl_count sqanti_abundance.csv \\
        --dir sqanti_out \\
        --output sqanti_output \\
        --cpus ${task.cpus} \\
        --report skip 

    # Ensure version check also uses the host python
    SQ_VERSION=\$(/env/SQANTI3/bin/python ${sqanti_dir}/sqanti3_qc.py --version 2>&1 | head -n1 | sed 's/SQANTI3 //')

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \${SQ_VERSION}
    END_VERSIONS
    """
}

process SQANTI3_REPORT {
    tag "SQANTI_Report"
    label 'process_single'
    
    input:
    path original_classification
    path junctions
    path sqanti_params

    output:
    path "*_report.html", emit: html
    path "*_report.pdf",  emit: pdf

    script:
    def sqanti_dir = "/data/gent/vo/000/gvo00082/vsc48905/Nextflow/apps/SQANTI3"
    def env_dir    = "/data/gent/vo/000/gvo00082/vsc48905/Nextflow/apps/sqanti_report_env"

    """
    # 1. The Ultimate HPC Conda Initializer (Bypasses .bashrc safety switches)
    if [ -f "\${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
        source "\${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "\${HOME}/miniforge3/etc/profile.d/conda.sh" ]; then
        source "\${HOME}/miniforge3/etc/profile.d/conda.sh"
    elif [ -f "\${VSC_DATA}/miniconda3/etc/profile.d/conda.sh" ]; then
        source "\${VSC_DATA}/miniconda3/etc/profile.d/conda.sh"
    else
        # Fallback to VSC module system if local installation isn't found
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load Miniconda3 2>/dev/null || module load Anaconda3 2>/dev/null || true
        eval "\$(conda shell.bash hook)" || true
    fi

    # 2. Build the working R environment (Only happens the first time!)
    if [ ! -d "${env_dir}" ]; then
        echo "Building persistent R environment for SQANTI3 reporting..."
        conda create -p "${env_dir}" -c conda-forge -c bioconda \
            r-base r-ggplot2 r-scales r-reshape r-reshape2 \
            r-dt r-knitr r-rmarkdown r-dplyr r-ggplotify \
            r-yulab.utils r-optparse r-gridextra bioconductor-noiseq r-plotly -y
    fi

    # 3. Safely run the host R script on the compute node
    conda run -p "${env_dir}" Rscript ${sqanti_dir}/utilities/report_qc/SQANTI3_report.R \\
        ${original_classification} \\
        ${junctions} \\
        ${sqanti_params} \\
        ${sqanti_dir}/utilities \\
        False \\
        both
    """
}