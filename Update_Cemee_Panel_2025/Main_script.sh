#!/bin/bash

### New script that uses Alix Desmarais's pipeline to include new sorted bam file to the Cemee panel ###
### It uses part of launch_GATK.part1.sh and launch_GATK.part2.sh found in her Github                ###
### It creates a new GVCF file with the input bam that can be then merged 


# Script which launch GATK steps until haplotype caller
# Path to change
path_ref='/projects/nematode_genetics/reference_genomes/20220216_c_elegans_WS283.genome.fa'
#path_gatk='/mnt/data2/mallard/2025_Cemee_panel_update/gatk-4.5.0.0/gatk'
path_gatk='/mnt/data2/desmarais/gatk-4.5.0.0/gatk'
# A clone of Alix git directory: git clone https://github.com/ExpEvolWormLab/Alix_internship.git
path_alix_pipeline='/mnt/data2/mallard/2025_Cemee_panel_update/Alix_internship/Pipeline_GATK/'

# location of the bam files to add
path_input_bam="/mnt/data2/mallard/2025_Cemee_panel_update/bam_sorted"
output_directory="/mnt/data2/mallard/2025_Cemee_panel_update/"

#path_known_sites='/mnt/data3/desmarais/Data/none.bed'
### Or create an empy bed file
touch /mnt/data2/mallard/2025_Cemee_panel_update/empty.bed
# index the empty bed file
$path_gatk IndexFeatureFile -I /mnt/data2/mallard/2025_Cemee_panel_update/empty.bed
path_known_sites='/mnt/data2/mallard/2025_Cemee_panel_update/empty.bed'

# Memory to allocate
#export JAVA_OPTS="-Xmx20g -XX:ConcGCThreads=5 -Djava.util.concurrent.ForkJoinPool.common.parallelism=10"

set -e
start=$(date +%s) #store time

#mkdir -p $output_directory
mkdir -p $output_directory/GATK

#### PREPROCESSING ####
# first feature : path of the file where .bam are stored : filtrerd_bam
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)

### Remove optical 

#chmod +x $path_alix_pipeline/script_preprocessing_f.sh
$path_alix_pipeline/script_preprocessing_f.sh $path_input_bam $output_directory/GATK $path_gatk $path_ref $path_known_sites 


#### VARIANT CALLING ####

# first feature : path of the file where .bam are stored : corrected_dedub.bam
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)
# third : chomosome number ( I II III IV V X) to call variant

path_2="$output_directory/GATK/BQSR_files"

# Array to hold the background process IDs
declare -a bg_pids
chmod +x $path_alix_pipeline/script_callingVariant_f.sh
for i in I II III IV V X MtDNA
do
	$path_alix_pipeline/script_callingVariant_f.sh $path_2 $output_directory/GATK $i $path_gatk $path_ref &
	# Store the background process ID
	bg_pids+=($!)
done

# Wait for all variant calling processes to finish
for pid in "${bg_pids[@]}"; do
    wait "$pid"
done

# Get the name of all sample
for i in $output_directory/GATK/HC_files/*.g.vcf.gz
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
        bcftools concat "$i.I.g.vcf.gz" "$i.II.g.vcf.gz" "$i.III.g.vcf.gz" "$i.IV.g.vcf.gz" "$i.V.g.vcf.gz" "$i.X.g.vcf.gz" "$i.MtDNA.g.vcf.gz" -O z -o $i.g.vcf.gz
	echo $i > sample.txt
	bcftools reheader -s sample.txt -O z -o $i.reheader.g.vcf.gz  $i.g.vcf.gz #Nem the sample with line ID
        $path_gatk IndexFeatureFile -I $i.reheader.g.vcf.gz
    else echo "PB with HC for $i"
    fi
done < list.name


#### END of the Launch_part1 from Alix

### Move the file to a new directory

mkdir reheader_files
mv GATK/HC_files/*reheader.g.vcf.gz* reheader_files/

##### START of the launch_part2  #####

### Create the .list file with all GVCF files

find -name GATK/HC_files/*reheader.g.vcf.gz > .list | sed "s/\.\//$pwd/g"

#### GENOTYPING ####
# first feature : path of the file where .list is stored
# second : name of the directory where the results'll be stored  (is going to be created where the script is launch)
#chmod +x $path_alix_pipeline/script_genotyping_f.1.sh

list_file="new_gvcf.list"
output_directory_genotyping="${output_directory}/genotyped_files"


### WARNING - Here add the path to the old GVCF file you want to update
ls reheader_files/*reheader.g.vcf.gz > $list_file
echo "/mnt/usb_FM/mallard/GVCFs/GVCF_Founders_AllLines.g.vcf.gz" >> $list_file
## Index the large GVCF if not done before
$path_gatk IndexFeatureFile -I "/mnt/usb_FM/mallard/GVCFs/GVCF_Founders_AllLines.g.vcf.gz"

## Local copy of the reference genome

cd $output_directory

### This create the gvcf file - all CHR at the same time, takes a few days to run
nohup  $path_alix_pipeline/script_genotyping_f.sh $list_file $output_directory_genotyping $path_gatk $path_ref "I" > ouput_genotype_chrI &
## For some reason it contains all CHR and then we should rename it /// Ideally modify the script to avoid this problem, simply remove .{chrom}
mv $output_directory_genotyping/All_combined.final.I.vcf.gz $output_directory_genotyping/All_combined.final.vcf.gz
mv $output_directory_genotyping/All_combined.final.I.vcf.gz.tbi $output_directory_genotyping/All_combined.final.vcf.gz.tbi

mv $output_directory_genotyping/All_combined.I.g.vcf.gz $output_directory_genotyping/All_combined.g.vcf.gz
mv $output_directory_genotyping/All_combined.I.g.vcf.gz.tbi $output_directory_genotyping/All_combined.g.vcf.gz.tbi


#### Then we go to SNP filtration #####

mkdir $output_directory/Filtered

# Keep only SNPs
## SelectVariant allow to only keep one nt as ref (for example if we got insertion as ATTA will only keep one)
$path_gatk SelectVariants -R $path_ref -V  $output_directory_genotyping/All_combined.final.vcf.gz --select-type-to-include SNP -O $output_directory/Filtered/combined.SNP.annotated.vcf.gz > SNP.filtration.output.txt

## Sometimes SelectVariant keep ref allele with only * : doesn't account as a SNP, remove them with :
bcftools view -v snps $output_directory/Filtered/combined.SNP.annotated.vcf.gz -O z > $output_directory/Filtered/combined.SNP.annotated_selected.vcf.gz
$path_gatk IndexFeatureFile -I $output_directory/Filtered/combined.SNP.annotated_selected.vcf.gz

# Using the threshold given by CaeNDR
$path_gatk VariantFiltration -R ${path_ref} -V $output_directory/Filtered/combined.SNP.annotated_selected.vcf.gz \
	--genotype-filter-expression "DP < 5.0"    --genotype-filter-name "DP_min_depth" \
            --filter-expression "QUAL < 30.0"                --filter-name "QUAL_quality" \
            --filter-expression "FS > 100.0"          --filter-name "FS_fisherstrand" \
            --filter-expression "QD < 20.0"      --filter-name "QD_quality_by_depth" \
            --filter-expression "SOR > 5.0"    --filter-name "SOR_strand_odds_ratio" \
            --genotype-filter-expression "isHet == 1.0"                  --genotype-filter-name "is_het" \
            -O $output_directory/Filtered/combined.annotated.filtered.vcf.gz > Variant.filtration.output.txt

# Apply high heterozygosity filter
bcftools filter --soft-filter='high_heterozygosity' --mode + --include 'F_PASS(GT="het")<=0.15' -O z  $output_directory/Filtered/combined.annotated.filtered.vcf.gz > $output_directory/Filtered/combined.annotated.filtered.het.vcf.gz
$path_gatk IndexFeatureFile -I $output_directory/Filtered/combined.annotated.filtered.het.vcf.gz

#Put ./. ie missing value for  heterozygotes genotype
$path_gatk SelectVariants \
	-V $output_directory/Filtered/combined.annotated.filtered.het.vcf.gz \
	--set-filtered-gt-to-nocall \
	-O $output_directory/Filtered/combined.without_hetero.vcf.gz > Variant.filtration.remove.het.genotypes.output.txt

# Apply high missing filter
bcftools filter --soft-filter='high_missing' --mode + --include 'F_MISSING<=0.95' $output_directory/Filtered/combined.without_hetero.vcf.gz -O z >  $output_directory/Filtered/combined.final.vcf.gz


### Write file which store SNP to exclude from imputation
# Multi-allelic SNP
bcftools view -H -m3 -v snps $output_directory/Filtered/combined.final.vcf.gz | awk '{print $1":"$2}' > $output_directory/Filtered/SNP_2exclude.txt
bcftools view -H -m3 -v snps $output_directory/Filtered/combined.final.vcf.gz | cut -f1,2 > $output_directory/Filtered/SNP_2exclude.tsv
echo "Multi-allelic SNP : "$(bcftools view -H -m3 -v snps $output_directory/Filtered/combined.final.vcf.gz |wc -l) > $output_directory/Filtered/summary.txt

# SNP with filter in non pass
bcftools view -H -f.,DP_min_depth,QUAL_quality,FS_fisherstrand,QD_quality_by_depth,SOR_strand_odds_ratio,high_missing,high_heterozygosity $output_directory/Filtered/combined.final.vcf.gz |\
 awk '{print $1":"$2}'  >> $output_directory/Filtered//SNP_2exclude.txt
bcftools view -H -f.,DP_min_depth,QUAL_quality,FS_fisherstrand,QD_quality_by_depth,SOR_strand_odds_ratio,high_missing,high_heterozygosity $output_directory/Filtered/combined.final.vcf.gz |\
 cut -f1,2  >> $output_directory/Filtered//SNP_2exclude.tsv
echo "non PASS SNP : "$(bcftools view -H -f.,DP_min_depth,QUAL_quality,FS_fisherstrand,QD_quality_by_depth,SOR_strand_odds_ratio,high_missing,high_heterozygosity $output_directory/Filtered/combined.final.vcf.gz | wc -l) >> $output_directory/Filtered/summary.txt


# Write singleton 
vcftools --gzvcf $output_directory/Filtered/combined.final.vcf.gz --singletons --out $output_directory/Filtered/to_remove

sed '1d' $output_directory/Filtered/to_remove.singletons | awk '{print $1":"$2}' >> $output_directory/Filtered/SNP_2exclude.txt
sed '1d' $output_directory/Filtered/to_remove.singletons | cut -f1,2 >> $output_directory/Filtered/SNP_2exclude.tsv
echo "doubletons : "$(sed '1d' $output_directory/Filtered/to_remove.singletons | wc -l) >> $output_directory/Filtered/summary.txt

####################################### Hard filtering ### Not necessary here - can be performed later, when unnecessary samples are removed/re-labeled #######################################
#bcftools view -T ^$output_directory/Filtered/SNP_2exclude.tsv $output_directory/Filtered/combined.final.vcf.gz -O z > $output_directory/Filtered/combined.final.hard_filter.vcf.gz
#########################################################################################################################################################################
######### Notes from Dehan LEE NEE paper  ###########################
#Finally, we removed sites that had >5% missing genotype data or where >10% of samples were called heterozygous.


