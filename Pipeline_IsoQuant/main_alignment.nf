// main.nf
nextflow.enable.dsl=2

include { LONGREAD_RNA } from './workflows/longread_rna_alignment'

workflow {
    LONGREAD_RNA()
}