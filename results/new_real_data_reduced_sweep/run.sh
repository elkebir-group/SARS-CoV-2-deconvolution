#/bin/bash

f="../../new_data_illumina_April8_reduced/final"

for k in 25 20 30 15 35 10 40 5 45 50
do
  for lambda in 1.005
  do
    for S in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
    do
      o=k${k}_S${S}
      echo "$1 -s $S -k $k -o $o -f ${f}_ref_results.tsv ${f}_alt_results.tsv --lambda ${lambda} --eps 0 --maxIter 200 > ${o}.log"
    done
  done
done
