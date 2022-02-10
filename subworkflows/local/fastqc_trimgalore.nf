//
// Read QC and trimming
//

params.fastqc_raw_options  = [:]
params.fastqc_trim_options = [:]
params.fastp_options       = [:]
trimgalore_options         = [:]

include { FASTQC as FASTQC_RAW  } from '../../modules/nf-core/modules/fastqc/main' addParams( options: params.fastqc_raw_options  )
include { FASTQC as FASTQC_TRIM } from '../../modules/nf-core/modules/fastqc/main' addParams( options: params.fastqc_trim_options )
include { TRIMGALORE       } from '../../modules/nf-core/modules/trimgalore/main'  addParams( options: params.trimgalore_options)


workflow FASTQC_TRIMGALORE {
    take:
    reads // channel: [ val(meta), [ reads ]
        
    main:
    fastqc_raw_html = Channel.empty()
    fastqc_raw_zip  = Channel.empty()
    fastqc_version  = Channel.empty()
    if (!params.skip_fastqc) {
        FASTQC_RAW ( reads ).html.set { fastqc_raw_html }
        fastqc_raw_zip = FASTQC_RAW.out.zip
        fastqc_version = FASTQC_RAW.out.versions
    }

    trim_reads       = reads
    trim_html        = Channel.empty()
    trim_log         = Channel.empty()
    trim_reads_fail  = Channel.empty()
    fastp_version    = Channel.empty()
    fastqc_trim_html = Channel.empty()
    fastqc_trim_zip  = Channel.empty()
    if (!params.skip_trimaglore) {
        TRIMGALORE (reads).reads.set { trim_reads }
        trim_log    = TRIMGALORE.out.log
        trimgalore_versions = TRIMGALORE.out.versions

        if (!params.skip_fastqc) {
            FASTQC_TRIM ( trim_reads ).html.set { fastqc_trim_html }
            fastqc_trim_zip = FASTQC_TRIM.out.zip
        }
    }

    emit:
    reads = trim_reads // channel: [ val(meta), [ reads ] ]
    trim_log           // channel: [ val(meta), [ log ] ]
   // trimgalore_version      //    path: *.version.txt

    fastqc_raw_html    // channel: [ val(meta), [ html ] ]
    fastqc_raw_zip     // channel: [ val(meta), [ zip ] ]
    fastqc_trim_html   // channel: [ val(meta), [ html ] ]
    fastqc_trim_zip    // channel: [ val(meta), [ zip ] ]
    fastqc_version     //    path: *.version.txt
}
