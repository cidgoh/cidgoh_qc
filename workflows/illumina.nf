/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

def modules = params.modules.clone()

def multiqc_options   = modules['illumina_multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

include { CUTADAPT              } from '../modules/local/cutadapt'              addParams( options: modules['illumina_cutadapt']      )
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files: ['tsv':'']]       )
include { MULTIQC               } from '../modules/local/multiqc_illumina'      addParams( options: multiqc_options                   )
include { KRAKEN2_DB_PREPARATION                              } from '../modules/local/kraken2_db_preparation'
include { KRAKEN2                                             } from '../modules/local/kraken2'                     addParams( options: modules['kraken2']                    )
include { KRONA_DB                                            } from '../modules/local/krona_db'
include { KRONA                                               } from '../modules/local/krona'                       addParams( options: modules['krona']                      )
include { CENTRIFUGE_DB_PREPARATION                           } from '../modules/local/centrifuge_db_preparation'
include { CENTRIFUGE                                          } from '../modules/local/centrifuge'  addParams( options: modules['centrifuge'] ) 
include { INPUT_CHECK        } from '../subworkflows/local/input_check'             addParams( options: [:] )


if(params.centrifuge_db){
    Channel
        .value(file( "${params.centrifuge_db}" ))
        .set { ch_centrifuge_db_file }
} else {
    ch_centrifuge_db_file = Channel.empty()
}

if(params.kraken2_db){
    Channel
        .value(file( "${params.kraken2_db}" ))
        .set { ch_kraken2_db_file }
} else {
    ch_kraken2_db_file = Channel.empty()
}


/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/
if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }
//
// MODULE: Installed directly from nf-core/modules
//
include { CAT_FASTQ                     } from '../modules/nf-core/modules/cat/fastq/main'                     addParams( options: modules['illumina_cat_fastq']                     )
include { FASTQC                        } from '../modules/nf-core/modules/fastqc/main'                        addParams( options: modules['illumina_cutadapt_fastqc']               )



//
// SUBWORKFLOW: Consisting entirely of nf-core/modules
//
def fastp_options = modules['illumina_fastp']
def trimgalore_options = modules['illumina_trimgalore']
if (params.save_trimmed_fail) { fastp_options.publish_files.put('fail.fastq.gz','seq/fail') }
if (params.save_trimmed_good) { fastp_options.publish_files.put('trim.fastq.gz','seq/trim') }
if (params.save_trimmed_good) { trimgalore_options.publish_files.put('fq.gz','seq/trim') }
if (params.save_merged) { fastp_options.publish_files.put('merged.fastq.gz','seq/trim') }

def bowtie2_align_options = modules['illumina_bowtie2_align']
if (params.save_unaligned) { bowtie2_align_options.publish_files.put('fastq.gz','unmapped') }

def markduplicates_options   = modules['illumina_picard_markduplicates']
markduplicates_options.args += params.filter_duplicates ?  Utils.joinModuleArgs(['REMOVE_DUPLICATES=true']) : ''


include { FASTQC_FASTP           } from '../subworkflows/local/fastqc_fastp'           addParams( fastqc_raw_options: modules['illumina_fastqc_raw'], fastqc_trim_options: modules['illumina_fastqc_trim'], fastp_options: fastp_options )
include { FASTQC_TRIMGALORE           } from '../subworkflows/local/fastqc_trimgalore'   addParams( fastqc_raw_options: modules['illumina_fastqc_raw'], fastqc_trim_options: modules['illumina_fastqc_trim'], trimgalore_options: trimgalore_options )


/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report    = []
def pass_mapped_reads = [:]
def fail_mapped_reads = [:]

workflow ILLUMINA {
    ch_versions = Channel.empty()

    INPUT_CHECK (
        ch_input,
        params.platform
    )
    .sample_info
    .map {
        meta, fastq ->
            meta.id = meta.id.split('_')[0..-2].join('_')
            [ meta, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)
    
    /*
    ================================================================================
             MODULE: Concatenate FastQ files from same sample if required
    ================================================================================
    */

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }


    if (params.qc_trim_mode == 'trimgalore' ) {
         FASTQC_TRIMGALORE (
          ch_cat_fastq
        )
        ch_kraken2_fastq    = FASTQC_TRIMGALORE.out.reads
        ch_centrifuge_fastq    = FASTQC_TRIMGALORE.out.reads

    }

    else {
        FASTQC_FASTP (
        ch_cat_fastq
        ) 
        ch_kraken2_fastq    =         FASTQC_FASTP.out.
        ch_centrifure_fastq    = FASTQC_TRIMGALORE.out.reads
        
        ch_versions = ch_versions.mix(FASTQC_FASTP.out.fastqc_version.first().ifEmpty(null))
        ch_versions = ch_versions.mix(FASTQC_FASTP.out.fastp_version.first().ifEmpty(null))
     
    }    

    /*
    ================================================================================
                                    Taxonomic information
    ================================================================================
    */
    CENTRIFUGE_DB_PREPARATION ( ch_centrifuge_db_file )
    CENTRIFUGE (
        ch_centrifuge_fastq,
        CENTRIFUGE_DB_PREPARATION.out.db
    )
    ch_versions = ch_versions.mix(CENTRIFUGE.out.version.first().ifEmpty(null))

    
    KRAKEN2_DB_PREPARATION (
        ch_kraken2_db_file
    )
    KRAKEN2 (
        ch_kraken2_fastq,
        KRAKEN2_DB_PREPARATION.out.db
        
    )
    ch_versions = ch_versions.mix(KRAKEN2.out.version.first().ifEmpty(null))
    
    if (( params.centrifuge_db || params.kraken2_db ) && !params.skip_krona) {
        KRONA_DB ()
        CENTRIFUGE.out.results_for_krona.mix(KRAKEN2.out.results_for_krona)
            . map { classifier, meta, report ->
                def meta_new = meta.clone()
                meta_new.classifier  = classifier
                [ meta_new, report ]
            }
            .set { ch_tax_classifications }
        KRONA (
            ch_tax_classifications,
            KRONA_DB.out.db.collect()
        )
        ch_versions = ch_versions.mix(KRONA.out.version.first().ifEmpty(null))
    }


    ch_versions.view()
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report, fail_mapped_reads)
    }
    NfcoreTemplate.summary(workflow, params, log, fail_mapped_reads, pass_mapped_reads)
}

/*
========================================================================================
    THE END
========================================================================================
*/
