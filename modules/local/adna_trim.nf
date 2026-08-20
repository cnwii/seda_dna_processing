process ADNA_TRIM {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::seqtk"

    container "/hpcfs/groups/acad_users/shyrav/resources/eager-v2.4.5-sharding-adna-trim-l25.sif"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.merged.fastq.gz")                      , emit: collapsed           , optional: true
    tuple val(meta), path("*_pe.R{1,2}.fq.gz")                      , emit: paired_truncated    , optional: true
    tuple val(meta), path('*.summary')                              , emit: summary

    script:
    def args = task.ext.args   ?: ''
    def base = "${reads[0].baseName}"
    """
    seqtk mergepe ${reads[0]} ${reads[1]} | \
    ${params.adna_trim_path} -t ${task.cpus} -p ${base}_pe ${args} - | gzip - > tmp.fq.gz

    seqkit sana tmp.fq.gz -o ${base}.merged.fastq.gz

    rm -f tmp.fq.gz
    cp .command.err ${base}_aDNA_trim.summary
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
seqtk mergepe $r1 $r2 | \
    ${params.adna_trim_path} -t ${task.cpus} -p ${base}_pe -l 25 - | gzip - > tmp.fq.gz

  seqkit sana tmp.fq.gz -o ${base}.merged.fastq.gz

  rm -f tmp.fq.gz

  cat ${base}.merged.fastq.gz ${base}_pe.R1.fq.gz ${base}_pe.R2.fq.gz > ${base}.adna_trim.fastq.gz

  cp .command.err ${base}_aDNA_trim.summary
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
