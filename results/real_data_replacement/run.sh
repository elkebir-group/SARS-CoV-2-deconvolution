#/bin/bash

for sample in *
do
	if [ -d $sample ]
	then
		#echo $sample
		sampleName=$sample
		for strain in $sampleName/*
		do
			#echo $strain 
			#exit 1
			if [ -d $strain ]
			then
				#echo $strain
				#exit 1
				refInput=$sampleName/ref_input.tsv
				altInput=$sampleName/alt_input.tsv
				
				for newBFile in $strain/*
				do
					caseName=$( basename $newBFile .tsv )
					if [ -f $newBFile ]
					then
						mkdir -p $strain/output
						outputPrefix=$strain/output/${caseName}_gradient
						#echo $refInput $altInput $newBFile $outputPrefix
						echo "../../build/gradient $refInput $altInput -B $newBFile -o $outputPrefix --lambda 0 -m 1 -k 25  > $outputPrefix.log"
						#exit 1
					fi
				done
			fi
		done
	fi
done
