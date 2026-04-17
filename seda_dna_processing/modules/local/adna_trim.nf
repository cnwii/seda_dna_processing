process ADNA_TRIM {
    tag "$meta.id"
    label 'process_medium'

    container "shyamrav/nf-core-eager:2.4.5-sharding"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${base}.merged.fastq.gz")                , emit: collapsed           , optional: true
    tuple val(meta), path("${base}_pe.R{1,2}.fq.gz")                , emit: paired_truncated    , optional: true
    tuple val(meta), path('*.summary')                              , emit: summary

    script:
    def args = task.ext.args   ?: ''
    if (meta.single_end) { 
        def base = "${reads.baseName}"
        """
        zcat $reads | ${params.adna_trim_path} -t ${task.cpus} ${args} - | gzip - > tmp.fq.gz
        seqkit sana tmp.fq.gz -o ${base}.merged.fastq.gz
        rm -f tmp.fq.gz
        cp .command.err ${base}_aDNA_trim.summary
        """
    } else {
        def base = "${reads[0].baseName}"
        """
        seqtk mergepe ${reads[0]} ${reads[1]} | \
        ${params.adna_trim_path} -t ${task.cpus} -p ${base}_pe ${args} - | gzip - > tmp.fq.gz

        seqkit sana tmp.fq.gz -o ${base}.merged.fastq.gz

        rm -f tmp.fq.gz
        cp .command.err ${base}_aDNA_trim.summary
        """
    }
}