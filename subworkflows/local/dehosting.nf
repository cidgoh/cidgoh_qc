
// import modules

include {generateCompositeReference} from '../modules/reads.nf'
include {grabCompositeIndex} from '../modules/reads.nf'
include {indexReference} from '../modules/reads.nf'
include {mapToCompositeIndex} from '../modules/reads.nf'
include {dehostBamFiles} from '../modules/reads.nf'
include {generateDehostedReads} from '../modules/reads.nf'
include {combineDehostedCSVs} from '../modules/reads.nf'




workflow vp_Dehosting {
    take:
      ch_filePairs
      ch_vpReference
      ch_HumanReference

    main:

      generateCompositeReference(ch_HumanReference, ch_vpReference)

      if ( params.bwa_index ){
        grabCompositeIndex("${params.bwa_index}")
        grabCompositeIndex.out
                .set{ ch_index }
      } else {
        indexReference(generateCompositeReference.out)
        indexReference.out
                .set{ ch_index }
      }

      mapToCompositeIndex(ch_filePairs.combine(generateCompositeReference.out.fasta),ch_index)

      dehostBamFiles(mapToCompositeIndex.out.bam)

      generateDehostedReads(dehostBamFiles.out.bam)

      generateDehostedReads.out.dehosted.set{ ch_filePairs_dehosted }

      combineDehostedCSVs(dehostBamFiles.out.csv.collect())

    emit:
      dehosted = ch_filePairs_dehosted

}
