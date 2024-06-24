#!/bin/bash

# Script which launch GATK steps until haplotype caller
# Path to change
path_ref='/mnt/data3/desmarais/Data/20220216_c_elegans_WS283.genome.fa'
path_gatk='/mnt/data3/desmarais/gatk-4.5.0.0/gatk'
# Memory to allocate
#export JAVA_OPTS="-Xmx30g -XX:ConcGCThreads=10 -Djava.util.concurrent.ForkJoinPool.common.parallelism=10"

set -e
start=$(date +%s) #store time

#### GENOTYPING ####
# first feature : path of the file where .list is stored
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)
for i in I II
do
	bash script_genotyping_f.1.sh $1 $2/GATK $path_gatk $path_ref $i &
done

#wait 

for i in III IV
do
        bash script_genotyping_f.1.sh $1 $2/GATK $path_gatk $path_ref $i &
done

wait 

for i in V X
do 
	bash script_genotyping_f.1.sh $1 $2/GATK $path_gatk $path_ref $i &
done

bash script_genotyping_f.1.sh $1 $2/GATK $path_gatk $path_ref MtDNA

wait 

# Concat all files together
bcftools concat $2/GATK/All_combined.final.*.vcf.gz -O z > $2/All_Lines_Chrom.final.vcf.gz

#### FILTRATION ####

# first feature : final file of previous HC
# second : name of the directory where the results'll be store
path_3="$2/All_Lines_Chrom.final.vcf.gz"
echo $path_3

bash script_filtration_f.sh $path_3 $2/GATK $path_gatk $path_ref


# Store end time
end=$(date +%s)
# Calculate running time
runtime=$((end - start))
# Echo running time
echo "launch_GATK.part2 executed successfully. Running time: $runtime seconds."
