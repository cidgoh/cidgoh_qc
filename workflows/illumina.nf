nextflow.enable.dsl=2



/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/


if (params.input)      { ch_input      = file(params.input)      } else { exit 1, 'Input samplesheet file not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

//ch_multiqc_config        = file("$projectDir/assets/multiqc_config_illumina.yaml", checkIfExists: true)
//ch_multiqc_config        = Channel.empty()
//ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

def modules = params.modules.clone()

//def multiqc_options   = modules['illumina_multiqc']
//multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/


include { INPUT_CHECK                 } from '../subworkflows/input_check'       addParams( options: [:] )
//include { GET_SOFTWARE_VERSIONS       } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
//include { MULTIQC                     } from '../modules/local/multiqc_illumina'      addParams( options: multiqc_options                   )
include { KRAKEN2_BUILD                     } from '../modules/local/kraken2_build'      addParams( options: [:]                   )




/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//


include { FASTQC                      } from '../modules/nf-core/modules/fastqc/main'    addParams([:])
include { TRIMGALORE                  } from '../modules/nf-core/modules/trimgalore/main'    addParams([:])
include { KRAKEN2_KRAKEN2             } from '../modules/nf-core/modules/kraken2/kraken2/main'    addParams([:])
include { FASTP                       } from '../modules/nf-core/modules/fastp/main'    addParams([:])
include { CUTADAPT                    } from '../modules/nf-core/modules/cutadapt/main'    addParams([:])
include { MULTIQC                     } from '../modules/nf-core/modules/multiqc/main'    addParams([:])

def multiqc_report    = []
workflow illumina{


    INPUT_CHECK (ch_input)



    FASTQC(INPUT_CHECK.out.sample_info)
    ch_fastqc_multiqc = FASTQC.out.zip



    ch_trimgalore_multiqc = Channel.empty()
    ch_cutadapt_multiqc = Channel.empty()
    ch_fastp_multiqc = Channel.empty()

    if (params.adapter_trim_mode == 'trimgalore' ) {
      ch_trimgalore    = INPUT_CHECK.out.sample_info
      TRIMGALORE(ch_trimgalore)
      ch_kraken2_fastq    = TRIMGALORE.out.reads
      ch_trimgalore_multiqc  = TRIMGALORE.out.log

    }
    else if (params.adapter_trim_mode == 'cutadapt' ){
      ch_cutadapt    = INPUT_CHECK.out.sample_info
      CUTADAPT(ch_cutadapt)
      ch_kraken2_fastq    = CUTADAPT.out.reads
      ch_cutadapt_multiqc  = CUTADAPT.out.log

    }

    else{
      save_trimmed_fail = false
      save_merged       = false
      ch_fastp    = INPUT_CHECK.out.sample_info
      FASTP(ch_fastp, save_trimmed_fail, save_merged)
      ch_kraken2_fastq    = FASTP.out.reads
      ch_fastp_multiqc  = FASTP.out.json

    }


    ch_kraken2_multiqc = Channel.empty()
    if (!params.skip_kraken2) {
        KRAKEN2_BUILD(params.kraken2_db)
        KRAKEN2_KRAKEN2 (
            ch_kraken2_fastq,
            KRAKEN2_BUILD.out.db
        )
        ch_kraken2_multiqc   = KRAKEN2_KRAKEN2.out.txt
    }

    if (!params.skip_multiqc) {
        //workflow_summary    = WorkflowCommons.paramsSummaryMultiqc(workflow, summary_params)
        //ch_workflow_summary = Channel.value(workflow_summary)

        MULTIQC (
            //ch_multiqc_config,
            //ch_multiqc_custom_config.collect().ifEmpty([]),
            //GET_SOFTWARE_VERSIONS.out.yaml.collect(),
            //ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),
            //ch_fail_reads_multiqc.ifEmpty([]),
            //ch_fail_mapping_multiqc.ifEmpty([]),
            //ch_amplicon_heatmap_multiqc.ifEmpty([]),
            ch_fastqc_multiqc.collect{it[1]}.ifEmpty([]),
            //ch_fastp_multiqc.collect{it[1]}.ifEmpty([]),
            //ch_kraken2_multiqc.collect{it[1]}.ifEmpty([]),
            //ch_cutadapt_multiqc.collect{it[1]}.ifEmpty([]),
        )
        multiqc_report = MULTIQC.out.report.toList()
    }


  }
