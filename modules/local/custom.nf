process generateCompositeReference {

    tag { "Vp_Human_CompositeIndex" }

    input:
    path(human_ref)
    path(vp_ref)

    output:
    path("composite_reference.fa"), emit: fasta

    script:
    """
    cat $human_ref $vp_ref > composite_reference.fa
    """
}

process grabCompositeIndex {

    tag { "Vp_Human_CompositeIndex" }

    input:
    path(index_folder)

    output:
    file("*.fa.*")

    script:
    """
    ln -sf $index_folder/*.fa.* ./
    """
}

process indexReference {

    tag { "bwa_composite_index" }

    input:
    path(composite_ref)

    output:
    file("*.fa*")

    script:
    """
    bwa index -a bwtsw $composite_ref
    """
}

process mapToCompositeIndex {
    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}.*", mode: "copy"

    //label 'mediumcpu'
    cpus 2

    input:
    tuple(sampleName, path(forward), path(reverse), path(composite_ref))
    path(indexed_reference)

    output:
    tuple sampleName, path("${sampleName}.sorted.bam"), emit: bam
    path("${sampleName}.flagstats.txt")

    script:
    """
    bwa mem -t ${task.cpus} ${composite_ref} ${forward} ${reverse} | samtools sort --threads ${task.cpus} -T "temp" -O BAM -o ${sampleName}.sorted.bam
    samtools flagstat ${sampleName}.sorted.bam > ${sampleName}.flagstats.txt
    """
}


process dehostBamFiles {
    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}*.{dehosted.bam,.csv}", mode: "copy"

    input:
    tuple(sampleName, path(composite_bam))

    output:
    tuple sampleName, path("${sampleName}.dehosted.bam"), emit: bam
    path("${sampleName}*.csv"), emit: csv

    script:

    def rev = workflow.commitId ?: workflow.revision ?: workflow.scriptId

    """
    samtools index ${composite_bam}
    dehost.py --file ${composite_bam} \
    -q ${params.keep_min_map_quality} \
    -Q ${params.remove_min_map_quality} \
    -o ${sampleName}.dehosted.bam \
    -R ${rev}
    """
}

process generateDehostedReads {

    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "${sampleName}-dehosted_R*", mode: "copy"

    input:
    tuple(sampleName, path(dehosted_bam))

    output:
    //tuple(sampleName, path("${sampleName}-dehosted_R1"), path("${sampleName}-dehosted_R2*")), emit: dehosted
    tuple sampleName, path("${sampleName}-dehosted_R1*"), path("${sampleName}-dehosted_R2*") , emit: dehosted

    script:
    """
    samtools fastq -1 ${sampleName}-dehosted_R1.fastq -2 ${sampleName}-dehosted_R2.fastq ${dehosted_bam}
    gzip ${sampleName}-dehosted_R*.fastq
    """
}

process combineDehostedCSVs {

    tag { "${params.prefix}/${sampleName}" }
    publishDir "${params.outdir}/${params.prefix}/${task.process.replaceAll(":","_")}", pattern: "*summary.csv", mode: "copy"

    input:
    path(csvs)

    output:
    path("removal_summary.csv")

    script:
    """
    csvtk concat *_stats.csv > removal_summary.csv
    """
}
