#/bin/bash

f="../../new_data_illumina_April8_reduced/final"

for k in 25 50
do
  for s in 0 1 2 3 4
  do
    o=mmf_k${k}_S${s}
    #grep 'normalized' ${o}.log | cut -d' ' -f 3 > ${o}.error 
    echo "$1 -s $s -k $k -o $o -f ${f}_ref_results.tsv ${f}_alt_results.tsv --lambda 0 --eps 0 --maxIter 200 > ${o}.log"
  done
done
