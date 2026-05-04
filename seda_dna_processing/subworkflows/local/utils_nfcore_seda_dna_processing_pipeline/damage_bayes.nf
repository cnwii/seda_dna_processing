process DAMAGE_BAYES {

    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.11 conda-forge::pandas conda-forge::numpy conda-forge::matplotlib conda-forge::pymc conda-forge::arviz"
    container "${ workflow.containerEngine == 'singularity' ?
    'docker://community.wave.seqera.io/library/biopython_matplotlib-base_numpy_pandas:fc3c9b3d64be2fab' :
    'community.wave.seqera.io/library/biopython_matplotlib-base_numpy_pandas:fc3c9b3d64be2fab' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.fastq'), emit: reads
    tuple val(meta), path('*.csv'), emit: damage_counts
    tuple val(meta), path('*.pdf'), emit: plot_damage

    script:
    """
    set -euo pipefail

    awk -v n_bp=${params.n_bp} -v trim=${params.trim} \
    -v out5="${reads.simpleName}_5_end_freq" \
    -v out3="${reads.simpleName}_3_end_freq" '
    BEGIN {
      print "Position_from_5end\\tA_freq\\tT_freq\\tC_freq\\tG_freq\\tTotal" > out5
      print "Position_from_3end\\tA_freq\\tT_freq\\tC_freq\\tG_freq\\tTotal" > out3
      }

    NR % 4 == 2 {
    seq = \$0
    len = length(seq)

    for (i = trim + 1; i <= n_bp; i++) {
      base5 = substr(seq, i, 1)
      base3 = substr(seq, len - i + 1, 1)

      count5[i,base5]++
      count3[i,base3]++

      total5[i]++
      total3[i]++
      }
    }

    END {
      for (i = trim + 1; i <= n_bp; i++) {
        printf "%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n", i, count5[i,"A"]+0, count5[i,"T"]+0, count5[i,"C"]+0, count5[i,"G"]+0, total5[i]+0 >> out5
        printf "%d\\t%d\\t%d\\t%d\\t%d\\t%d\\n", i, count3[i,"A"]+0, count3[i,"T"]+0, count3[i,"C"]+0, count3[i,"G"]+0, total3[i]+0 >> out3
      }
    }
    ' "${reads}"

    mkdir -p base_frequency plot_damage
    mv *fastq.gz base_frequency
    damage_bayes.py base_frequency/ plot_damage/
    """
}