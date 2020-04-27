#/bin/bash

for sample in *
do
	if [ -d $sample ]
	then
		#echo $sample
		sampleName=$sample
		refInput=$sampleName/ref_input.tsv
		altInput=$sampleName/alt_input.tsv
		Bfile='Bexpansion_input.tsv'	
		for k in $( seq 1 10 )
		do
			echo "../../build/exposure -k $k $refInput $altInput -B $Bfile -o ${sampleName}/exposure_$k > ${sampleName}/exposure_${k}.log"
		done
	fi
done
