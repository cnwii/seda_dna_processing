process KRAKEN_FILTERING {

    tag "${meta.id}_k${param_set.uniq_kmer}_r${param_set.tax_reads}"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::numpy"
    container "/hpcfs/groups/acad_users/containers/ngspy_0.3.sif"

    input:
    tuple val(meta), path(report), val(param_set)

    output:
    tuple val(meta), val(param_set), path("formatted/*formatted.report.txt.filtered"), emit: formatted

    script:
    def uniq_kmer = param_set.uniq_kmer
    def tax_reads = param_set.tax_reads

    def only_ratio_flag = params.only_ratio ? "--only_ratio" : ""
    def only_reads_flag = params.only_reads ? "--only_reads" : ""

    """
    set -euo pipefail

    new_report=\$(basename "${report}" | sed "s/\\.report/.k${param_set.uniq_kmer}_r${param_set.tax_reads}.formatted.report/")

    awk 'BEGIN{FS=OFS="\\t"}
    {
        gsub(/\\r/, "")
        sub(/[ \\t]+\$/, "")
        print
    }' "${report}" > "\$new_report"

    krakenuniq_filter.py \\
        --krakenuniq_report \$new_report \\
        --n_unique_kmers ${uniq_kmer} \\
        --n_tax_reads ${tax_reads} \\
        --ratio ${params.ratio} \\
        --rank ${params.rank} \\
        ${only_ratio_flag} ${only_reads_flag}

    mkdir -p formatted

    mv *formatted.report.txt.filtered formatted/
    """
}

process KRAKEN_ABUNDANCES {

    label 'process_single'

    conda "${projectDir}/subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/filtering_environment.yml"

    input:
    val(meta)
    tuple val(meta2), val(param_set), path(formatted_file)

    output:
    tuple val(meta), val(param_set), path("abundance/*.txt"), emit: abundance

    script:
    """
    set -euo pipefail

    mkdir -p formatted abundance

    cp ${formatted_file} formatted

    krakenuniq_abundances.R \\
    formatted \\
    abundance \\
    ${meta.id} \\
    ${param_set.uniq_kmer} \\
    ${param_set.tax_reads}
    """
}

process KRAKEN_PLOT {

    label 'process_single'

    conda "${projectDir}/subworkflows/local/utils_nfcore_seda_dna_processing_pipeline/filtering_environment.yml"

    input:
    tuple val(meta), val(param_set), path(abundance_files)
    path(tax_file)
    path(metadata)

    output:
    tuple val(meta), val(param_set), val("taxonomy"), path("plots/*"), emit: plots
    tuple val(meta), val(param_set), path("diversity/*"), emit: diversity


    script:
    """
    set -euo pipefail

    mkdir -p abundance plots diversity
    cp ${abundance_files} abundance

    taxonomy_profiles.R \\
    abundance \\
    plots \\
    ${tax_file} \\
    ${meta.id} \\
    ${metadata} \\
    ${param_set.uniq_kmer} \\
    ${param_set.tax_reads}

    beta_diversity.R \\
    abundance \\
    diversity \\
    ${tax_file} \\
    ${meta.id} \\
    ${metadata} \\
    ${param_set.uniq_kmer} \\
    ${param_set.tax_reads} 
    """
}