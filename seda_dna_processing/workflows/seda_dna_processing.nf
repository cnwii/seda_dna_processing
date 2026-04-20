
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK                                   } from '../subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/input_check'
include { FASTQC                                        } from '../modules/nf-core/fastqc/main'
include { FASTP as FASTP_POLYG                          } from '../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PAIRED       } from '../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SINGLE       } from '../modules/nf-core/adapterremoval/main'
include { ADNA_TRIM                                     } from '../modules/local/adna_trim.nf'
include { CAT_FASTQ as CAT_FASTQ_AR                     } from '../modules/nf-core/cat/fastq/main'
include { FASTP as FASTP_LOW_COMPLEXITY                 } from '../modules/nf-core/fastp/main'
//include { BBMAP_BBDUK                                   } from '../modules/nf-core/bbmap/bbduk/main' 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEDA_DNA_PROCESSING {

    take:
    samplesheet // channel: samplesheet read in from --input
    
    main:

    // Check input
    INPUT_CHECK (samplesheet)

    // Run basic QC
    FASTQC (
        INPUT_CHECK.out.reads
    )

    // Trim PolyG tail
    FASTP_POLYG (
        INPUT_CHECK.out.reads.map { meta, reads ->
        tuple(meta, reads, [])
    },
    false, false, false
    )

    // Branch PE and SE data
    ch_input_for_adapterremoval = FASTP_POLYG.out.reads
        .branch{
            single: it[0]['single_end'] == true
            paired: it[0]['single_end'] == false
        }
    // Run AR for SE
    ADAPTERREMOVAL_SINGLE (
        ch_input_for_adapterremoval.single, file(params.adapter_list)
    )

    if(!params.use_adna_trim) {
        // Run AR for PE
        ADAPTERREMOVAL_PAIRED (
            ch_input_for_adapterremoval.paired, file(params.adapter_list)
        )

        // Take collapsed and collapsed Truncated FASTQs from AR_PE, and concatinate them into on FASTQ
        // Convert the PE reads which have been collapsed to SE. 
        ch_concat_ar_pe_fastq = ADAPTERREMOVAL_PAIRED.out.collapsed
            .mix(ADAPTERREMOVAL_PAIRED.out.collapsed_truncated)
            .map {
                meta, reads -> 
                    def meta_new = meta.clone()
                    meta_new.single_end = true
                    [ meta_new, reads ]
            }
            .groupTuple()
            .map { meta, reads -> [ meta, reads.flatten() ] }

        CAT_FASTQ_AR(ch_concat_ar_pe_fastq)

        // Take only Collapsed PE and Truncated SE for further processing
        ch_adapter_trimmed_reads_prepped = CAT_FASTQ_AR.out.reads
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

    } else {
        // ADNA_TRIM 
        ADNA_TRIM(ch_input_for_adapterremoval.paired)

        // Convert the PE reads which have been collapsed to SE. 
        ch_adapter_trimmed_reads_prepped = ADNA_TRIM.out.collapsed
            .map {
                meta, reads -> 
                    def meta_new = meta.clone()
                    meta_new.single_end = true
                    [ meta_new, reads ]
            }
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)
    }

    
    FASTP_LOW_COMPLEXITY (
        ch_adapter_trimmed_reads_prepped.map { meta, reads ->
        tuple(meta, reads, [])
        },
        false, false, false
    )

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/