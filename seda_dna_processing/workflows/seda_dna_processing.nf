
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK                                   } from '../subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/input_check'
include { FASTQC as FASTQC_RAW_READS                    } from '../modules/nf-core/fastqc/main'
include { FASTP as FASTP_POLYG                          } from '../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PAIRED       } from '../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_SINGLE       } from '../modules/nf-core/adapterremoval/main'
include { ADNA_TRIM                                     } from '../modules/local/adna_trim.nf'
include { CAT_FASTQ as CAT_FASTQ_AR                     } from '../modules/nf-core/cat/fastq/main'
include { FASTP as FASTP_LOW_COMPLEXITY                 } from '../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_PROCESSED_READS              } from '../modules/nf-core/fastqc/main'
include { UNTAR as UNTAR_METAGENOMICS                   } from '../modules/nf-core/untar/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ                } from '../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'
include { KRAKENTOOLS_KREPORT2KRONA                     } from '../modules/nf-core/krakentools/kreport2krona/main'                 

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEDA_DNA_PROCESSING {

    take:
    samplesheet // channel: samplesheet read in from --input
    
    main:

// Pre-processing
    // Check input
    INPUT_CHECK (samplesheet)

    // Run basic QC
    FASTQC_RAW_READS (
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

    // Filtering low complexity reads
    FASTP_LOW_COMPLEXITY (
        ch_adapter_trimmed_reads_prepped.map { meta, reads ->
        tuple(meta, reads, [])
        },
        false, false, false
    )

    // QC report for processed reads
    FASTQC_PROCESSED_READS (
        FASTP_LOW_COMPLEXITY.out.reads
    )

    // Metagenomics profiling
    ch_database = Channel.value(file(params.database_path))
    .branch{
        untar: it ==~ /.*.tar.gz/
        base:true
    }

    // Untar the database
    ch_untar_input = ch_database.untar.map{ [[], it] }

    UNTAR_METAGENOMICS ( 
        ch_untar_input 
        )

    ch_untar_output = UNTAR_METAGENOMICS.out.untar.map{ it[1] }

    ch_database = ch_database.base.mix(ch_untar_output)

    ch_krakenuniq_input = FASTP_LOW_COMPLEXITY.out.reads
    .map { meta, reads ->
        def seqtype = meta.single_end ? 'fastq' : 'fastq'
        def prefix  = meta.id

        tuple(
            [meta, reads, [prefix]],
            seqtype,
            null                      
        )
    }
    .combine(ch_database)
    .map { tuple_part, seqtype, _, db ->
        def (meta, reads, prefixes) = tuple_part
        tuple(
            [meta, reads, prefixes],
            seqtype,
            db
        )
    }

    // Run KrakenUniq
    KRAKENUNIQ_PRELOADEDKRAKENUNIQ(
        ch_krakenuniq_input.map{ it[0] },
        ch_krakenuniq_input.map{ it[1] },
        ch_krakenuniq_input.map{ it[2] },
        true,
        true,
        true
    )

    KRAKENTOOLS_KREPORT2KRONA(
        KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/