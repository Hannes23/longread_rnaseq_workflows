// workflows/longread_rna.nf

// Importing modules
include { FASTQC } from '../modules/local/fastqc'
include { FASTPLONG } from '../modules/local/fastplong'
include { MINIMAP2 } from '../modules/local/minimap2'
include { CLEAN_GENOME } from '../modules/local/transcriptclean'
include { SPLIT_BAM_BY_CHR } from '../modules/local/split_bam'
include { SAMTOOLS_MERGE } from '../modules/local/samtools_merge'
include { TRANSCRIPT_CLEAN } from '../modules/local/transcriptclean'
include { TALON_INIT; TALON_RUN; TALON_CREATE_GTF } from '../modules/local/talon'
include { SQANTI3; SQANTI3_REPORT } from '../modules/local/sqanti3'
include { SQANTI3_FILTER; SQANTI3_FILTER_REPORT } from '../modules/local/sqanti3_filter'
include { SALMON_INDEX; SALMON_QUANT } from '../modules/local/salmon'
include { MULTIQC } from '../modules/local/multiqc'

// Workflow definition
workflow LONGREAD_RNA {

    // 1. VALIDATE INPUT
    if (params.input == null) {
        exit 1, "ERROR: Please provide a samplesheet via --input samplesheet.csv"
    }

    // 2. PARSE SAMPLESHEET
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id = row.sample
            def fastq = file(row.fastq)
            if (!fastq.exists()) {
                exit 1, "ERROR: FastQ file not found: ${row.fastq}"
            }
            return [ meta, fastq ]
        }
        .set { ch_reads }

    // 3. REFERENCES & ASSETS
    ch_genome = Channel.value(file(params.genome))
    ch_gtf    = Channel.value(file(params.gtf))
    
    // NEW: Define the path to your strict JSON filter rules
    ch_filter_rules = Channel.value(file("${projectDir}/assets/filtering.json"))

    // --- PIPELINE LOGIC ---
    
    // Clean genome
    CLEAN_GENOME(ch_genome)

    // QC
    FASTQC(ch_reads)
    FASTPLONG(ch_reads)

    // Alignment
    MINIMAP2(FASTPLONG.out.reads, ch_genome)

    // Split by chromosome
    SPLIT_BAM_BY_CHR(MINIMAP2.out.bam)

    // Flatten output to (meta, bam) tuples
    ch_chrom_bams = SPLIT_BAM_BY_CHR.out.chromosome_bams
        .flatMap { meta, bam_list ->
            bam_list.collect { bam ->
                // Parse chromosome name from filename (e.g. sample_chr1.bam -> chr1)
                // Assuming format: "${meta.id}_${chrom}.bam"
                def chrom_name = bam.name.toString()
                                    .minus("${meta.id}_")
                                    .minus(".bam")
                
                // Return tuple matching TranscriptClean input
                tuple(meta, chrom_name, bam)
            }
        }

    // Transcript cleaning per chromosome
    // Use .first() to ensure genome is ready and available
    TRANSCRIPT_CLEAN(
        ch_chrom_bams,
        CLEAN_GENOME.out.clean_genome.first() 
    )
    
    // Merge cleaned BAMs by sample
    ch_cleaned_per_sample = TRANSCRIPT_CLEAN.out.bam
        .groupTuple(by: 0)
        .map { meta, bam_list ->
            tuple(meta, bam_list)
        }

    SAMTOOLS_MERGE(ch_cleaned_per_sample)

    // TALON discovery
    TALON_INIT(ch_gtf)

    // Flatten merged BAM channel to just file paths for TALON_RUN
    ch_merged_bams = SAMTOOLS_MERGE.out.bam
        .map { tuple -> tuple[1] }  // Only BAM path
        .collect()

    TALON_RUN(
        ch_merged_bams,
        TALON_INIT.out.db,
        ch_genome, 
        params.platform
    )
    
    // TALON EXPORT
    TALON_CREATE_GTF(
        TALON_RUN.out.db
    )

    // SQANTI3 QC
    SQANTI3(
        TALON_CREATE_GTF.out.gtf,       
        TALON_CREATE_GTF.out.abundance,
        ch_gtf,
        ch_genome
    )

    SQANTI3_REPORT(
        SQANTI3.out.original_classification,
        SQANTI3.out.junctions,
       SQANTI3.out.sqanti_params
    )

    // SQANTI3 FILTER (UPDATED)
    // Now accepts the JSON rules file as the 4th argument
    SQANTI3_FILTER(
        SQANTI3.out.original_classification,
        SQANTI3.out.fasta,
        SQANTI3.out.corrected_gtf,
        ch_filter_rules   // <--- Passed the JSON file here
    )

      SQANTI3_FILTER_REPORT(
        SQANTI3.out.original_classification,        
        SQANTI3_FILTER.out.filtered_classification, 
        SQANTI3_FILTER.out.reasons         
    )

    // Salmon quantification
    SALMON_INDEX(
        ch_genome,                  
        SQANTI3_FILTER.out.fasta          
    )

  
    SALMON_QUANT(
        FASTPLONG.out.reads,
        SALMON_INDEX.out.index
    )

    // MultiQC
     // Gather all QC-relevant outputs
    ch_multiqc_files = Channel.empty()

    // FASTQC
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip)
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.html)

    // FASTPLONG
    ch_multiqc_files = ch_multiqc_files.mix(FASTPLONG.out.html)

    // SQANTI3 reports
    ch_multiqc_files = ch_multiqc_files.mix(SQANTI3_REPORT.out.report)

    // SQANTI3 filter reports
    ch_multiqc_files = ch_multiqc_files.mix(SQANTI3_FILTER_REPORT.out.report)

    // SALMON outputs (quant.sf + logs)
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.quant)
    ch_multiqc_files = ch_multiqc_files.mix(SALMON_QUANT.out.logs)

    // Run MultiQC once on all collected files
    MULTIQC(
        ch_multiqc_files.collect()
    )

}