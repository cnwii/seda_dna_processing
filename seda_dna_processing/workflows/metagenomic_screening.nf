/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BWA_INDEX         } from '../modules/nf-core/bwa/index/main' 
include { BWA_ALN           } from '../modules/nf-core/bwa/aln/main' 
include { BOWTIE2_BUILD     } from '../modules/nf-core/bowtie2/build/main' 
include { BOWTIE2_ALIGN     } from '../modules/nf-core/bowtie2/align/main'
include { BLAST_BLASTN      } from '../modules/nf-core/blast/blastn/main'

include { MALT_BUILD        } from '../modules/nf-core/malt/build/main'   
include { MALT_RUN          } from '../modules/nf-core/malt/run/main' 
include { MALTEXTRACT       } from '../modules/nf-core/maltextract/main'
include { MEGAN_RMA2INFO    } from '../modules/nf-core/megan/rma2info/main'  
include { MAPDAMAGE2        } from '../modules/nf-core/mapdamage2/main'                                                                                       

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAGENOMIC_SCREENING {

    take:
    samplesheet 
    
    main:

    

    BWA_INDEX(

    )

    BWA_ALN(

    )

    BLAST_BLASTN(

    )



}