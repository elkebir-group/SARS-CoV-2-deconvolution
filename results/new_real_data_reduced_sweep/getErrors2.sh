#/bin/bash


for k in 25 50
do
  for s in 0 1 2 3 4
  do
    o=mmf_k${k}_S${s}
    grep 'normalized' ${o}.log | cut -d' ' -f 3 > ${o}.error
  done
done
