                                                                                                                                                    #!/usr/bin/env nextflow
/*
========================================================================================
    cidgoh_qc
========================================================================================
    Github : https://github.com/cidgoh/cidgoh_qc
    Website: https://cidgoh.ca/
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    GENOME PARAMETER VALUES
========================================================================================
*/


params.fasta         = WorkflowMainQC.getGenomeAttribute(params, 'fasta')


/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/

WorkflowMainQC.initialise(workflow, params, log)

/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

if (params.platform == 'illumina') {
    include { ILLUMINA } from './workflows/illumina'
} else if (params.platform == 'nanopore') {
    include { NANOPORE } from './workflows/nanopore'
}


//
// WORKFLOW: Run main QC analysis pipeline
//

workflow CIDGOH_QC {
    //
    // WORKFLOW: QC analysis for Illumina data
    //
    if (params.platform == 'illumina') {
        ILLUMINA ()

    //
    // WORKFLOW: QC analysis for Nanopore data
    //
    } else if (params.platform == 'nanopore') {
        NANOPORE ()
    }
}


/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/


workflow {
    CIDGOH_QC ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
