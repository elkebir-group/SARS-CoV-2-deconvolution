#/bin/bash

echo sample k error
for sample in *
do
	if [ -d $sample ]
	then
		for k in $( seq 1 10 )
		do
			error=$( grep 'Frob norm' ${sample}/exposure_${k}.log | cut -d' ' -f 7 )
			echo $sample $k $error
		done
	fi
done
