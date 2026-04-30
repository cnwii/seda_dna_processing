process KRAKEN_FILTERING {

    tag "${meta.id}_k${param_set.uniq_kmer}_r${param_set.tax_reads}"
    label 'process_single'

    conda """
    conda-forge::python=3.8.3
    conda-forge::pandas
    """
    container "/hpcfs/groups/acad_users/containers/ngspy_0.3.sif"

    input:
    tuple val(meta), path(report), val(param_set)

    output:
    tuple val(meta), val(param_set), path("*.txt"), emit: formatted_report

    script:
    def uniq_kmer = param_set.uniq_kmer
    def tax_reads = param_set.tax_reads

    def only_ratio_flag = params.only_ratio ? "--only_ratio" : ""
    def only_reads_flag = params.only_reads ? "--only_reads" : ""

    """
    set -euo pipefail

    new_report=\$(basename ${report} | sed 's/.report/.k${param_set.uniq_kmer}_r${param_set.tax_reads}.formatted.report/')

    sed -e "s/\\r//g" ${report} > \$new_report

    krakenuniq_filter.py \\
        --krakenuniq_report \$new_report \\
        --n_unique_kmers ${uniq_kmer} \\
        --n_tax_reads ${tax_reads} \\
        --ratio ${params.ratio} \\
        --rank ${params.rank} \\
        ${only_ratio_flag} ${only_reads_flag}
    """
}

process KRAKEN_PLOT {

    tag "${meta.id}_k${param_set.uniq_kmer}_r${param_set.tax_reads}"
    label 'process_single'

    conda "conda-forge::r-base=4.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/r-base:4.3.1' :
        'biocontainers/r-base:4.3.1' }"

    input:
    tuple val(meta), val(param_set), path(formatted_report)

    output:
    tuple val(meta), val(param_set), path("formatted/*"), emit: formatted
    tuple val(meta), val(param_set), path("abundance/*"), emit: abundance
    tuple val(meta), val(param_set), path("plots/*"),     emit: plots

    script:
    """
    set -euo pipefail

    mkdir -p formatted abundance plots

    cp -f -- ${formatted_report} formatted/

    krakenuniq_abundances.R formatted abundance
    plot_abundances.R abundance plots
    """
}