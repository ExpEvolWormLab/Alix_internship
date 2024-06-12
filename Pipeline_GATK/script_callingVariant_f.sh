#!/bin/bash

# first feature : path of the file where .bam are stored
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)
# third : chomosome number ( I II III IV V X) to call variant

path_gatk=$4
path_ref=$5


set -e #Stop if error
start=$(date +%s) #store time

# Create directories
mkdir -p $2
mkdir -p $2/HC_files

for file in "$1"/*.bam
do
	line=${file%_*}
	line=${line##*/}

	# Calling variants
	for chrom in $3
	do
	$path_gatk HaplotypeCaller --java-options "-Xmx20g" \
            --emit-ref-confidence GVCF \
            -R ${path_ref} \
            -I ${file} \
		-L $chrom \
	--annotation DepthPerAlleleBySample \
               --annotation Coverage \
               --annotation GenotypeSummaries \
               --annotation TandemRepeat \
           --annotation StrandBiasBySample \
          --annotation ChromosomeCounts \
            --annotation ReadPosRankSumTest \
            --annotation AS_ReadPosRankSumTest \
            --annotation AS_QualByDepth \
               --annotation QualByDepth \
            --annotation AS_StrandOddsRatio \
            --annotation AS_MappingQualityRankSumTest \
            --annotation DepthPerSampleHC \
            --annotation-group StandardAnnotation \
            --annotation-group AS_StandardAnnotation \
            --annotation-group StandardHCAnnotation \
            -O $2/HC_files/${line}.${chrom}.g.vcf.gz \
		 --tmp-dir $2/tmp

	done
done

# Store end time
end=$(date +%s)
# Calculate running time
runtime=$((end - start))
# Echo running time
echo "Script executed successfully. Running time: $runtime seconds."
