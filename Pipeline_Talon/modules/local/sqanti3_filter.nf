process SQANTI3_FILTER {

    tag "$original_classification"
    label 'process_medium'

    input:
    path original_classification
    path corrected_fasta
    path corrected_gtf
    path filter_json

    output:
    path "sqanti_filter_out/*.filtered.gtf",                emit: gtf
    path "sqanti_filter_out/*.filtered.fasta",              emit: fasta
    path "sqanti_filter_out/*classification.txt",           emit: filtered_classification
    path "sqanti_filter_out/*reasons.txt",                  emit: reasons
    path "versions.yml",                                    emit: versions

    script:
    def sqanti_dir = "/data/gent/vo/000/gvo00082/vsc48905/Nextflow/apps/SQANTI3"
    def env_dir    = "/data/gent/vo/000/gvo00082/vsc48905/Nextflow/apps/sqanti_report_env"

    """
    set -euo pipefail

    # 1. Initialize Conda for the compute node
    if [ -f "\${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
        source "\${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "\${HOME}/miniforge3/etc/profile.d/conda.sh" ]; then
        source "\${HOME}/miniforge3/etc/profile.d/conda.sh"
    elif [ -f "\${VSC_DATA}/miniconda3/etc/profile.d/conda.sh" ]; then
        source "\${VSC_DATA}/miniconda3/etc/profile.d/conda.sh"
    else
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load Miniconda3 2>/dev/null || module load Anaconda3 2>/dev/null || true
        eval "\$(conda shell.bash hook)" || true
    fi

    # 2. Build environment ONLY if it doesn't already exist
    if [ ! -d "${env_dir}" ]; then
        conda create -p "${env_dir}" -c conda-forge -c bioconda \
            python=3.11 r-base r-ggplot2 r-scales r-reshape r-reshape2 \
            r-dt r-knitr r-rmarkdown r-dplyr r-ggplotify \
            r-yulab.utils r-optparse r-gridextra bioconductor-noiseq r-plotly pandas pysam bx-python biopython -y
    fi

    # THE FIX: Added the base cDNA_Cupcake directory to the PYTHONPATH!
    export PYTHONPATH="${sqanti_dir}/cDNA_Cupcake/sequence/:${sqanti_dir}/cDNA_Cupcake/:\${PYTHONPATH:-}"

    mkdir -p sqanti_filter_out

    echo "1. Adding Roy score..."
    conda run -p "${env_dir}" python \\
        /data/gent/vo/000/gvo00082/vsc48905/Nextflow/Pipeline/bin/add_roy_score.py \\
        ${original_classification} \\
        scored_classification.txt

    echo "2. Running SQANTI3 rules filter..."
    conda run -p "${env_dir}" python \\
        ${sqanti_dir}/sqanti3_filter.py rules \\
        scored_classification.txt \\
        -j ${filter_json} \\
        --isoforms ${corrected_fasta} \\
        --gtf ${corrected_gtf} \\
        -d sqanti_filter_out \\
        -o sqanti_filtered \\
        --skip_report

    echo "3. Capturing version..."
    # Safe YAML generation block
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: 5.2
    END_VERSIONS
    """
}


process SQANTI3_FILTER_REPORT {

    tag "SQANTI_Filter_Report"
    label 'process_single'

    input:
    path original_classification
    path filtered_classification
    path reasons

    output:
    path "*report.html",         emit: html, optional: true
    path "*report.pdf",          emit: pdf
    path "versions.yml",         emit: versions

    script:
    def sqanti_dir = "/data/gent/vo/000/gvo00082/vsc48905/Nextflow/apps/SQANTI3"
    def env_dir    = "/data/gent/vo/000/gvo00082/vsc48905/Nextflow/apps/sqanti_report_env"

    """
    set -euo pipefail

    # 1. Initialize Conda for the compute node
    if [ -f "\${HOME}/miniconda3/etc/profile.d/conda.sh" ]; then
        source "\${HOME}/miniconda3/etc/profile.d/conda.sh"
    elif [ -f "\${HOME}/miniforge3/etc/profile.d/conda.sh" ]; then
        source "\${HOME}/miniforge3/etc/profile.d/conda.sh"
    elif [ -f "\${VSC_DATA}/miniconda3/etc/profile.d/conda.sh" ]; then
        source "\${VSC_DATA}/miniconda3/etc/profile.d/conda.sh"
    else
        source /etc/profile.d/modules.sh 2>/dev/null || true
        module load Miniconda3 2>/dev/null || module load Anaconda3 2>/dev/null || true
        eval "\$(conda shell.bash hook)" || true
    fi

    # 2. Stage the filter results properly
    mkdir -p filter_results
    
    # We copy the files into the directory to preserve the EXACT names SQANTI3 generated
    cp ${filtered_classification} filter_results/
    cp ${reasons} filter_results/

    echo "Running SQANTI3 filter report..."

    # FIXED: Removed the invalid '-c' flag!
    conda run -p "${env_dir}" Rscript \\
        ${sqanti_dir}/utilities/report_filter/SQANTI3_filter_report.R \\
        -d filter_results \\
        -o sqanti_filtered \\
        -u ${sqanti_dir}/utilities \\
        -f rules

    echo "Capturing R version..."
    R_VERSION=\$(conda run -p "${env_dir}" Rscript --version 2>&1 | head -n1)

    # Safe YAML generation block
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \${R_VERSION}
    END_VERSIONS
    """
}