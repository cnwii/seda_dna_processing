
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
include { DAMAGE_BAYES                                  } from '../subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/damage_bayes'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ                } from '../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'
include { KRAKEN_FILTERING                              } from '../subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/krakenuniq_filtering'
include { KRAKEN_ABUNDANCES                             } from '../subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/krakenuniq_filtering'
include { KRAKEN_PLOT                                   } from '../subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/krakenuniq_filtering'
include { KRAKENTOOLS_KREPORT2KRONA                     } from '../modules/nf-core/krakentools/kreport2krona/main'
include { MULTIQC as MULTIQC_RAW                        } from '../modules/nf-core/multiqc/main'
include { MULTIQC as MULTIQC_PROCESSED                  } from '../modules/nf-core/multiqc/main'
//include { QUARTONOTEBOOK                                } from '../modules/nf-core/quartonotebook/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SEDA_DNA_PROCESSING {

    take:
    samplesheet // channel: samplesheet read in from --input
    
    main:

// PRE-PROCESSING
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

    //MULTIQC runs     
    ch_raw_reads = Channel.empty()
    ch_processed_reads = Channel.empty()

    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)

    ch_raw_reads = FASTQC_RAW_READS.out.zip
    .map { meta, file -> file }
    .collect()

    ch_processed_reads = FASTQC_PROCESSED_READS.out.zip
    .map { meta, file -> file }
    .collect()

    ch_multiqc_meta = FASTQC_PROCESSED_READS.out.zip
    .map { meta, file -> meta }


    MULTIQC_RAW(
        ch_multiqc_meta,
        ch_raw_reads,
        ch_multiqc_config,
        [],
        [],
        []
    )

    MULTIQC_PROCESSED(
        ch_multiqc_meta,
        ch_processed_reads,
        ch_multiqc_config,
        [],
        [],
        []
    )

// ESTIMATING DAMAGE (DAMAGE BAYES)
    DAMAGE_BAYES (
        FASTP_LOW_COMPLEXITY.out.reads
    )


// METAGENOMICS PROFILING
    ch_database = Channel.fromPath(params.database_path, type: 'dir')

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

    // KrakenUniq filtering
    report_ch = KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report
// Set multiple thresholds
    //
    param_combos = Channel.of(
        [uniq_kmer: 1000, tax_reads: 100],
        [uniq_kmer: 500,  tax_reads: 50]
    )

    // attach params to each upstream output
    report_ch
        .combine(param_combos)
        .map { meta, report, param_set ->
            tuple(meta, report, param_set)
            }
        .set { input_ch }


    KRAKEN_FILTERING(input_ch)
    
    ch_tax_file = Channel.fromPath(params.tax_file)
    ch_metadata = Channel.fromPath(params.metadata)

    ch_id = KRAKEN_FILTERING.out.formatted
    .map { meta, param_set, file ->
        meta
    }
    ch_reports = KRAKEN_FILTERING.out.formatted
        .groupTuple(by:1)

    
    KRAKENTOOLS_KREPORT2KRONA(
        KRAKEN_FILTERING.out.formatted  
        )

    KRAKEN_ABUNDANCES(
        ch_id,
        ch_reports
    )

    ch_abundances = ch_abundances = KRAKEN_ABUNDANCES.out.abundance
    .groupTuple(by:1)
    .map { metas, param_set, abundance_files ->
        tuple(
            metas[0],
            param_set,
            abundance_files.flatten()
        )
    }
    ch_abundances.view()


    KRAKEN_PLOT(
        ch_abundances,
        ch_tax_file,
        ch_metadata
    )

    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GRAVEYARD

    ch_quarto_input = ch_kraken_meta
    .combine(ch_quartonotebook)
    .map { meta_plots, notebook ->
        def (meta, plots) = meta_plots
        tuple(meta, notebook)
    } // Give meta to notebook

    KRAKEN_PLOT.out.plots
    .map { meta, param_set, taxonomy, plots ->
        tuple(taxonomy, plots)
    }
    .mix()

    // Create an empty .qmd notebook beforehand for input
    
    ch_quartonotebook = Channel.fromPath(params.quartonotebook)

    ch_quarto_input = KRAKEN_FILTERING.out.formatted
    .combine(ch_quartonotebook)
    .map { meta, param_set, file, notebook ->
        tuple(meta, notebook)
    }

    ch_params_block = [
        dna_concentration: "number",
        library_fragments: "number",
        concentration: "value",
        analysis: "method"
    ]

    ch_input_pdf = 
        DAMAGE_BAYES.out.plot_damage
            .map { meta, plots ->
                plots
            }
            .collect()

    
    ch_kraken_meta = KRAKEN_PLOT.out.plots.map { meta, param_set, taxonomy, plots ->
    tuple(meta, plots)
    }

QUARTONOTEBOOK(
        ch_quarto_input,
        ch_params_block,
        ch_input_pdf,
        []
        )

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/