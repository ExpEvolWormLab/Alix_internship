#!/bin/bash

# first feature : path of the folder where .bam are stored
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)

path_gatk=$3
path_ref=$4
path_known_sites=$5


set -e
start=$(date +%s)

# Create directories
mkdir -p $2
# Store the sort files
mkdir -p $2/sort_files
# Store the flagged duplicates files
mkdir -p $2/dedub_files
# Store base quality score recalibration table and results
mkdir -p $2/BQSR_files

mkdir -p $2/tmp

for file in "$1"/*.bam
do
        #Get the name of the line
        line=${file%%.*}
        line=${line##*/}
	if grep -q ${line} $(ls $2/BQSR_files/*);
	then
		continue
	fi
	echo -e "LINE : $line\n\n"
 
        #Flag the duplicates
        $path_gatk MarkDuplicates --java-options "-Xmx20g" \
				-I $file \
                                -O $2'/dedub_files/'$line'.sorted.dedub.bam' \
                                -M $2'/dedub_files/'$line'.sorted.dedub.metrics.txt' \
				--TMP_DIR $2/tmp
    
        #BQSR
	## Add ReadGroup if the file doesn't have any
	echo -e "PATH : $2'/dedub_files/'$line'.sorted.dedub.bam'\n\n\n"
	echo "FILE : $file"
	echo -e "\n\n\n\n"
	if ! samtools view -H $file | grep -q '@RG'; then
		$path_gatk AddOrReplaceReadGroups --java-options "-Xmx20g" \
				        -I $2'/dedub_files/'$line'.sorted.dedub.bam' \
					-O $2'/dedub_files/'$line'.sorted.dedub.RG.bam' \
					-RGID $line \
				        -RGLB lib1 \
       					-RGPL ILLUMINA \
			         	-RGPU unit1 \
				        -RGSM $line \
					--TMP_DIR $2/tmp
	else
		cp $2'/dedub_files/'$line'.sorted.dedub.bam' $2'/dedub_files/'$line'.sorted.dedub.RG.bam'
	fi

        ## Compute recalibaration table
	samtools index -m 8G $2'/dedub_files/'$line'.sorted.dedub.RG.bam'
        $path_gatk BaseRecalibrator --java-options "-Xmx20g" \
				-I $2'/dedub_files/'$line'.sorted.dedub.RG.bam'\
                                -R $path_ref \
                                --known-sites $path_known_sites \
                                -O $2'/BQSR_files/'$line'_recal_data.table' \
				 --tmp-dir $2/tmp
        ## Correct BQS
        $path_gatk ApplyBQSR --java-options "-Xmx20g" \
				-I $2'/dedub_files/'$line'.sorted.dedub.RG.bam'\
                                -R $path_ref \
                                --bqsr-recal-file $2'/BQSR_files/'$line'_recal_data.table' \
                                -O $2'/BQSR_files/'$line'.sorted.dedub.corrected.bam' \
				 --tmp-dir $2/tmp
	rm $2'/dedub_files/'$line'.sorted.dedub.RG.bam'
	rm $2'/dedub_files/'$line'.sorted.dedub.RG.bam.csi'
	rm $2'/dedub_files/'$line'.sorted.dedub.bam'
	rm $2'/dedub_files/'$line'.sorted.dedub.metrics.txt'
	rm $2'/BQSR_files/'$line'_recal_data.table'
done
# Store end time
end=$(date +%s)
# Calculate running time
runtime=$((end - start))
# Echo running time
echo "Script executed successfully. Running time: $runtime seconds."


