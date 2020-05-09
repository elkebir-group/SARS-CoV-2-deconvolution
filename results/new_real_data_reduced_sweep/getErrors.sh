#/bin/bash

for k in 5 10 15 20 25 30 35 40 45 50
do
  for lambda in 1.005
  do
    for S in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
    do
      o=k${k}_S${S}
      grep 'normalized' ${o}.log | cut -d' ' -f 3 > ${o}.error
    done
  done
done
