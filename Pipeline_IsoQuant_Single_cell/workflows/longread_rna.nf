// Importing modules
include { FASTQC } from '../modules/local/fastqc'
include { FLEXIPLEX_DISCOVER; FLEXIPLEX_FILTER; FLEXIPLEX_DEMUX } from '../modules/local/flexiplex'
include { MINIMAP2 } from '../modules/local/minimap2'
include { TAG_BAM } from '../modules/local/tag_bam'
include { SAMTOOLS_MERGE } from '../modules/local/samtools_merge'
include { ISOQUANT } from '../modules/local/isoquant'
include { SQANTI3 } from '../modules/local/sqanti3'
include { SQANTI3_FILTER } from '../modules/local/sqanti3_filter'

// Workflow definition
workflow LONGREAD_RNA {

    // 1. VALIDATE INPUT
    if (params.input == null) {
        exit 1, "ERROR: Please provide a samplesheet via --input samplesheet.csv"
    }

    // --- DEFINE REFERENCE CHANNELS --- 
    ch_genome       = Channel.value(file(params.genome))
    ch_gtf          = Channel.value(file(params.gtf))
    data_type       = params.data_type
    ch_filter_rules = Channel.value(file("${projectDir}/assets/filtering.json"))

    // 2. PARSE SAMPLESHEET
    ch_base_reads = Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .map { row ->
            def meta = [:]
            meta.id    = row.sample
            meta.group = row.condition
            def fastq = file(row.fastq)
            return [ meta, fastq ]
        }
    
    // 3. CHUNK THE FASTQ FILES
    ch_fastq_chunks = ch_base_reads.splitFastq(by: 10000000, file: true)

    // --- PIPELINE ---
    // Run FastQC on the WHOLE file
    FASTQC( ch_base_reads )

    // --- Flexiplex ---
    FLEXIPLEX_DISCOVER(ch_base_reads, params.sc_chemistry)
    FLEXIPLEX_FILTER(FLEXIPLEX_DISCOVER.out.counts)

    ch_demux_input = ch_fastq_chunks.combine(FLEXIPLEX_FILTER.out.final_barcodes, by: 0)

    FLEXIPLEX_DEMUX(
        ch_demux_input,
        params.sc_chemistry,
        params.sc_platform
    )

    // --- Alignment ---
    MINIMAP2(
        FLEXIPLEX_DEMUX.out.demux_fastq,
        ch_genome,
        data_type
    )

    // Inject the single-cell tags
    TAG_BAM(
        MINIMAP2.out.bam
    )

    // --- Prepare IsoQuant Input ---
    // Group BAMs back together per sample (e.g. sample1_chunk1, sample1_chunk2)
    ch_bams_to_merge = TAG_BAM.out.bam
        .map { meta, bam, bai -> tuple(meta, bam) }   
        .groupTuple(by: 0) 

    SAMTOOLS_MERGE( ch_bams_to_merge )


    ch_isoquant_input = SAMTOOLS_MERGE.out.merged_bam
        .map { meta, bam, bai -> 
            tuple(meta.group, bam, bai, meta.id)
        }
        .groupTuple() 

    // --- Quantification ---
    ISOQUANT(
        ch_isoquant_input,
        ch_genome,
        ch_gtf,
        params.data_type,
        params.umi_tag,
        params.barcode_tag,
        params.sc_mode
    )

    // --- SQANTI3  ---
    SQANTI3(
        ISOQUANT.out.gtf,
        ISOQUANT.out.counts,
        ch_gtf,
        ch_genome
    )

    // --- FILTER ---
    ch_classification = SQANTI3.out.original_classification
        .map { file ->
            def group_name = file.baseName.replace('_sqanti_classification', '')
            tuple(group_name, file)
        }

    SQANTI3_FILTER(
        ch_classification,
        SQANTI3.out.fasta,
        SQANTI3.out.corrected_gtf,
        ch_filter_rules
    )

}