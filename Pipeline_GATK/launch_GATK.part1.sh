#!/bin/bash

# Script which launch GATK steps until haplotype caller
# Path to change
path_ref='/mnt/data3/desmarais/Data/20220216_c_elegans_WS283.genome.fa'
path_known_sites='/mnt/data3/desmarais/Data/none.bed'
path_gatk='/mnt/data3/desmarais/gatk-4.5.0.0/gatk'
# Memory to allocate
#export JAVA_OPTS="-Xmx20g -XX:ConcGCThreads=5 -Djava.util.concurrent.ForkJoinPool.common.parallelism=10"


set -e
start=$(date +%s) #store time

mkdir -p $2
mkdir -p $2/GATK

#### PREPROCESSING ####
# first feature : path of the file where .bam are stored : filtrerd_bam
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)

bash script_preprocessing_f.sh $1 $2/GATK $path_gatk $path_ref $path_known_sites 


#### VARAINT CALLING ####

# first feature : path of the file where .bam are stored : corrected_dedub.bam
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)
# third : chomosome number ( I II III IV V X) to call variant

path_2="$2/GATK/BQSR_files"

# Array to hold the background process IDs
declare -a bg_pids

#for i in I II III IV V X MtDNA
for i in MtDNA
do
	bash script_callingVariant_f.sh $path_2 $1/GATK $i $path_gatk $path_ref &
	# Store the background process ID
	bg_pids+=($!)
done

# Wait for all variant calling processes to finish
for pid in "${bg_pids[@]}"; do
    wait "$pid"
done

# Get the name of all sample
for i in *.g.vcf.gz
	do echo ${i%%.*} >> list.name 
done

cat list.name | uniq > list.name1

rm list.name
mv list.name1 list.name 

# Rename each .gz with line as sample name
# Parcourir chaque ligne du fichier list.name
while IFS= read -r i; do
    # Vérifier si les fichiers existent pour chaque modèle de nommage
    if [ -e "$i.I.g.vcf.gz" ] && [ -e "$i.II.g.vcf.gz" ] && [ -e "$i.III.g.vcf.gz" ] && [ -e "$i.IV.g.vcf.gz" ] && [ -e "$i.V.g.vcf.gz" ] && [ -e "$i.X.g.vcf.gz" ] && [ -e "$i.MtDNA.g.vcf.gz" ]; then
        bcftools concat "$i.I.g.vcf.gz" "$i.II.g.vcf.gz" "$i.III.g.vcf.gz" "$i.IV.g.vcf.gz" "$i.V.g.vcf.gz" "$i.X.g.vcf.gz" "$i.MtDNA.g.vcf.gz "-O z -o $i.g.vcf.gz
	echo $i > sample.txt
	bcftools reheader -s sample.txt -O z -o $i.reheader.g.vcf.gz  $i.g.vcf.gz #Nem the sample with line ID
        $gatk IndexFeatureFile -I $i.reheader.g.vcf.gz
    else echo "PB with HC for $i"
    fi
done < list.name


# Store end time
end=$(date +%s)
# Calculate running time
runtime=$((end - start))
# Echo running time
echo "launch_GATK.part1 executed successfully. Running time: $runtime seconds."
echo "Need to create .list before running launch_GATK.part2."
echo 'find -name *reheader.g.vcf.gz > .list | sed "s/\.\//$pwd/g"'
