#/bin/bash
for k in 5 10
do
  for M in 0 0.1
  do
    for s in {0..20}
    do
      f=../../data/simulations/s${s}_M${M}_k${k}
      for S in {0..20}
      do
        o=s${s}_M${M}_k${k}_S${S}
        echo "$1 -s $S -k $k -o $o ${f}_ref.tsv ${f}_alt.tsv -l 1.005 --eps 0 --maxIter 100 > ${o}.log"
      done
    done
  done
done
