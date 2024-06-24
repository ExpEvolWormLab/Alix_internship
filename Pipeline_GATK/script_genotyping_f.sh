#!/bin/bash

# first feature : path of the file where .table is stored
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)


path_gatk=$3
path_ref=$4
chrom=$5


set -e #Stop if error
start=$(date +%s) #store time

mkdir -p $2
mkdir -p $2/tmp

#Merged all the vcf files
$path_gatk  CombineGVCFs\
       -R ${path_ref} \
        --variant $1 \
        -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
        -O $2/All_combined.${chrom}.g.vcf.gz \
	 --tmp-dir $2/tmp

#Perform joint genotyping on samples 

$path_gatk --java-options "-Xmx10g" GenotypeGVCFs \
   -R ${path_ref} \
   -V $2/All_combined.${chrom}.g.vcf.gz \
   -O $2/All_combined.final.${chrom}.vcf.gz \
	 --tmp-dir $2/tmp

# Store end time
end=$(date +%s)
# Calculate running time
runtime=$((end - start))
# Echo running time
echo "Script executed successfully. Running time: $runtime seconds."
