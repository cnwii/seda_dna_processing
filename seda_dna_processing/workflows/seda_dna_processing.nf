
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK            } from '../subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/input_check'

include { FASTQC                 } from '../modules/nf-core/fastqc/main'
include { FASTP                  } from '../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL         } from '../modules/nf-core/adapterremoval/main'
include { BBMAP_BBDUK            } from '../modules/nf-core/bbmap/bbduk/main'  

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEDA_DNA_PROCESSING {

    take:
    samplesheet // channel: samplesheet read in from --input
    
    main:

    INPUT_CHECK (samplesheet)


    FASTQC (
        INPUT_CHECK.out.reads
    )

    
    FASTP (
        INPUT_CHECK.out.reads.map { meta, reads ->
        tuple(meta, reads, [])
    },
    false, false, false
    )

    //FASTP.out.reads.view()

    adapter_list = Channel.value(file('/hpcfs/users/a1844642/project/seda_dna_processing/data/adapters.txt'))

    ADAPTERREMOVAL (
        FASTP.out.reads, adapter_list
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/