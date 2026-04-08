#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC } from '../modules/nf-core/fastqc/main'
include { FASTP } from '../modules/nf-core/fastp/main'
include { ADAPTERREMOVAL } from '../modules/nf-core/adapterremoval/main'
// include { ADNA_TRIM } from '../modeules/local/adna_trim/main'

include { BBMAP_BBDUK } from '../modules/nf-core/bbmap/bbduk/main' // This is low complexity reads
// include { BBMAP_DEDUPE } from '../modeules/local/bbmap/dedupe/main' // This is dedudping reads

include { VSEARCH_DEREPLICATE } from '../modules/nf-core/vsearch/dereplicate/main'

// include { DEDUP } from '../modules/nf-core/dedup/main' // This is deduping bams.

// include { KRAKEN2_KRAKEN2 } from '../modules/nf-core/kraken2/kraken2/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK } from '../subworkflow/input_check.nf'

workflow TEST {

    ch_samplesheet = Channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true)


    ch_samplesheet.view()

        INPUT_CHECK (
        file(ch_samplesheet)
    )

    FASTQC (
        INPUT_CHECK.out.reads
    )

    FASTP (
        INPUT_CHECK.out.reads, false, false, false
    )

    //ADAPTERREMOVAL(FASTP.out.reads, false)

    //BBMAP_BBDUK(ADAPTERREMOVAL.out.reads)
    //VSEARCH_DEREPLICATE(BBMAP_BBDUK.out.reads)
    
}

workflow {
    TEST()
}

//if (param.run_adnatrim) {
//        ADNA_TRIM(ch_input)
//        ch_collapsed_reads = ADNA_TRIM.out
//    } else {
//        ADAPTERREMOVAL(ch_input)
//        ch_collapsed_reads = ADAPTERREMOVAL.out
//    }