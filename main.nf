#!/usr/bin/env nextflow

/* 
 * pipeline input parameters 
 */

params.input = "$baseDir/sample_data/*_R{1,2}.fastq.gz"
params.outdir = "${PWD}/results"

def helpMessage() {
    log.info qcHeader()
    log.info """\
         CIDGOH - QC PIPELINE (v0.1)  
         ===================================
         input        : ${params.input}
         outdir       : ${params.outdir}
         """         .stripIndent()
}


def qcHeader(){
    return """
        ░█████╗░██╗██████╗░░██████╗░░█████╗░██╗░░██╗
        ██╔══██╗██║██╔══██╗██╔════╝░██╔══██╗██║░░██║
        ██║░░╚═╝██║██║░░██║██║░░██╗░██║░░██║███████║
        ██║░░██╗██║██║░░██║██║░░╚██╗██║░░██║██╔══██║
        ╚█████╔╝██║██████╔╝╚██████╔╝╚█████╔╝██║░░██║
        ░╚════╝░╚═╝╚═════╝░░╚═════╝░░╚════╝░╚═╝░░╚═╝
    """.stripIndent()
}

helpMessage()

Channel 
    .fromFilePairs(params.input, checkIfExists: true )
    .set{ read_pairs_ch }

process fastqc {
    tag "FASTQC"
    publishDir "${params.outdir}/fastqc_report"

    input:
    tuple sample_id, path(reads) from read_pairs_ch

    output:
    path "fastqc_${sample_id}_logs" into fastqc_ch

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  

process multiqc {
    tag "MULTIQC"
    publishDir params.outdir, mode:'copy'
       
    input:
    path '*' from fastqc_ch.collect()
    
    output:
    path 'multiqc_report.html'
     
    script:
    """
    multiqc . 
    """
} 



workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Script finished SUCCESSFULLY after $workflow.duration"
    } else {
      log.info "[$workflow.complete] >> Script finished with ERRORS after $workflow.duration"
    }
}