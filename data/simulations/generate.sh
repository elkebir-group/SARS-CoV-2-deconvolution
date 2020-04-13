#/bin/bash
for k in 5 10
do
  for M in 0 0.1
  do
    for s in {0..20}
    do
      $1 -s $s --missing $M -k $k -o s${s}_M${M}_k${k}
    done
  done
done
