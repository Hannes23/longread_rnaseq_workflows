// workflows/longread_rna.nf

// Importing modules
include { FASTQC } from '../modules/local/fastqc'
include { FASTPLONG } from '../modules/local/fastplong'
include { MINIMAP2 } from '../modules/local/minimap2'
include { ISOQUANT } from '../modules/local/isoquant'
include { SQANTI3 } from '../modules/local/sqanti3'
include { SQANTI3_FILTER } from '../modules/local/sqanti3_filter'


include { PREPARE_ISOQUANT_REFERENCE } from '../subworkflows/local/prepare_isoquant_reference'
// Workflow definition
workflow LONGREAD_RNA {

    // 1. VALIDATE INPUT
    if (params.input == null) {
        exit 1, "ERROR: Please provide a samplesheet via --input samplesheet.csv"
    }

    // 2. PARSE SAMPLESHEET
    Channel
        .fromPath(params.input)
        .splitCsv(header: true, strip: true)
        .map { row ->

            def meta = [:]
            meta.id        = row.sample.trim()
            meta.group = row.condition.trim()

            def fastq = file(row.fastq, checkIfExists: true)

            return [ meta, fastq ]
        }
        .set { ch_reads }
    // 3. REFERENCES & ASSETS
    ch_genome = Channel.value(file(params.genome))
    ch_gtf    = Channel.value(file(params.gtf))
    data_type = params.data_type
    ch_filter_rules = Channel.value(file("${projectDir}/assets/filtering.json"))

    // --- PIPELINE LOGIC ---

    // 1. QC
    FASTQC(ch_reads)
    FASTPLONG(ch_reads)

    // 2. Alignment
    MINIMAP2(FASTPLONG.out.reads, ch_genome, data_type)

    // 3. Directly feed each BAM per sample into IsoQuant
    ch_samples = MINIMAP2.out.bam
        .map { meta, bam, bai -> tuple(meta.group, bam, bai, meta.id) }
        .set { ch_sample_inputs }

    PREPARE_ISOQUANT_REFERENCE(
        ch_gtf
    )
    

    // 4. IsoQuant 
    ISOQUANT(
        ch_sample_inputs,
        ch_genome,
        PREPARE_ISOQUANT_REFERENCE.out.db.first(), // <--- THE CRITICAL FIX
        data_type
    )
    
    
    ch_sqanti_gtf = ISOQUANT.out.gtf.map { group, gtf_files ->

        def files = gtf_files instanceof List ? gtf_files : [gtf_files]

        def extended = files.find {
            it.name.endsWith("extended_annotation.gtf")
        }

        tuple(group, extended)
    }

    ch_sqanti_abundance = ISOQUANT.out.counts.map { group, count_files ->

        def files = count_files instanceof List ? count_files : [count_files]

        def abundance = files.find {
            it.name.endsWith("transcript_tpm.tsv")
        }

        tuple(group, abundance)
    }

    ch_for_sqanti = ch_sqanti_gtf.join(ch_sqanti_abundance)
    
    // 5. SQANTI3
    SQANTI3(
        ch_for_sqanti,
        ch_gtf,
        ch_genome
    )

    ch_classification = SQANTI3.out.original_classification
    .map { meta, file ->
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
