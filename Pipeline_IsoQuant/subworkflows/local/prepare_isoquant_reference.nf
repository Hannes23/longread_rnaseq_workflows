include { ISOQUANT_DB } from '../../modules/local/isoquant_db'

workflow PREPARE_ISOQUANT_REFERENCE {
    take:
    gtf

    main:
    // If the user provided a DB file via command line, use it directly
    if (params.provided_db) {
        db = Channel.fromPath(params.provided_db)
    } 
    // Otherwise, run the process to create it
    else {
        db = ISOQUANT_DB(gtf)
    }

    emit:
    db
}