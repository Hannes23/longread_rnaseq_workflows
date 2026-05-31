// main.nf
nextflow.enable.dsl=2

include { LONGREAD_RNA } from './workflows/longread_rna'

workflow {
    LONGREAD_RNA()
}