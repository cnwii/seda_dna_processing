#!/bin/bash

# Calculate nucleotide frequencies
# Parameters
N_bp=25
TRIM=1

Count_ATCG() {
  FILE=$1

  awk -v n_bp="$N_bp" -v trim="$TRIM" \
      -v out5="${FILE/.fastq}_5_end_freq" \
      -v out3="${FILE/.fastq}_3_end_freq" '
    BEGIN {
      print "Position_from_5end\tA_freq\tT_freq\tC_freq\tG_freq\tTotal" > out5
      print "Position_from_3end\tA_freq\tT_freq\tC_freq\tG_freq\tTotal" > out3
    }

    NR % 4 == 2 {
      seq = $0
      len = length(seq)

      for (i = trim + 1; i <= n_bp; i++) {
        base5 = substr(seq, i, 1)
        base3 = substr(seq, len - i + 1, 1)

        count5[i][base5]++
        count3[i][base3]++

        total5[i]++
        total3[i]++
      }
    }

    END {
      for (i = trim + 1; i <= n_bp; i++) {
        printf "%d\t%d\t%d\t%d\t%d\t%d\n", i, count5[i]["A"]+0, count5[i]["T"]+0, count5[i]["C"]+0, count5[i]["G"]+0, total5[i]+0 >> out5
        printf "%d\t%d\t%d\t%d\t%d\t%d\n", i, count3[i]["A"]+0, count3[i]["T"]+0, count3[i]["C"]+0, count3[i]["G"]+0, total3[i]+0 >> out3
      }
    }
  ' "$FILE"
}

for FILE in *.fastq; do
  echo "Processing $FILE..."
  Count_ATCG "$FILE"
done

echo "Base frequency calculation complete."