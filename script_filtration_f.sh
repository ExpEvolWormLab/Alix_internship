#!/bin/bash

# first feature : final file of previous HC
# second : name of the directory where the results'll be store

# To change :
path_ref=$4
path_gatk=$3



set -e #Stop if error
start=$(date +%s) #store time

# Create directories
mkdir -p $2/Filtred

# Variant Filtration
# Keep only SNPs
## SelectVariant allow to only keep one nt as ref (for example if we got insertion as ATTA will only keep one)
$path_gatk SelectVariants \
     -R $path_ref \
     -V  $1 \
     --select-type-to-include SNP \
     -O $2/Filtred/combined.SNP.annotated.vcf.gz

## Sometimes SelectVariant keep ref allele with only * : doesn't account as a SNP, remove them with :
bcftools view -v snps $2/Filtred/combined.SNP.annotated.vcf.gz -O z > $2/Filtred/combined.SNP.annotated.1.vcf.gz
$path_gatk IndexFeatureFile -I $2/Filtred/combined.SNP.annotated.1.vcf.gz

# I'm using the threshold given by CaeNDR : not optimal
$path_gatk VariantFiltration \
   -R ${path_ref}\
   -V $2/Filtred/combined.SNP.annotated.1.vcf.gz \
	--genotype-filter-expression "DP < 5.0"    --genotype-filter-name "DP_min_depth" \
            --filter-expression "QUAL < 30.0"                --filter-name "QUAL_quality" \
            --filter-expression "FS > 100.0"          --filter-name "FS_fisherstrand" \
            --filter-expression "QD < 20.0"      --filter-name "QD_quality_by_depth" \
            --filter-expression "SOR > 5.0"    --filter-name "SOR_strand_odds_ratio" \
            --genotype-filter-expression "isHet == 1.0"                  --genotype-filter-name "is_het" \
            -O $2/Filtred/combined.annotated.filtred.vcf.gz


# Apply high heterozygosity filter
bcftools filter --soft-filter='high_heterozygosity' --mode + --include 'F_PASS(GT="het")<=0.15' -O z  $2/Filtred/combined.annotated.filtred.vcf.gz > $2/Filtred/combined.annotated.filtred.1.vcf.gz
$path_gatk IndexFeatureFile -I $2/Filtred/combined.annotated.filtred.1.vcf.gz

#Put ./. ie missing value for  heterozygotes genotype
$path_gatk SelectVariants \
	-V $2/Filtred/combined.annotated.filtred.1.vcf.gz \
	--set-filtered-gt-to-nocall \
	-O $2/Filtred/combined.without_hetero.vcf.gz

# Apply high missing filter
bcftools filter --soft-filter='high_missing' --mode + --include 'F_MISSING<=0.95' $2/Filtred/combined.without_hetero.vcf.gz -O z >  $2/Filtred/combined.final.vcf.gz


### Write file which store SNP to exclude from imputation
# Multi-allelic SNP
bcftools view -H -m3 -v snps $2/Filtred/combined.final.vcf.gz | awk '{print $1":"$2}' > $2/Filtred/SNP_2exclude.txt
bcftools view -H -m3 -v snps $2/Filtred/combined.final.vcf.gz | cut -f1,2 > $2/Filtred/SNP_2exclude.tsv
echo "Multi-allelic SNP : "$(bcftools view -H -m3 -v snps $2/Filtred/combined.final.vcf.gz |wc -l) > $2/Filtred/summary.txt

# SNP with filter in non pass
bcftools view -H -f.,DP_min_depth,QUAL_quality,FS_fisherstrand,QD_quality_by_depth,SOR_strand_odds_ratio,high_missing,high_heterozygosity $2/Filtred/combined.final.vcf.gz |\
 awk '{print $1":"$2}'  >> $2/Filtred/SNP_2exclude.txt
bcftools view -H -f.,DP_min_depth,QUAL_quality,FS_fisherstrand,QD_quality_by_depth,SOR_strand_odds_ratio,high_missing,high_heterozygosity $2/Filtred/combined.final.vcf.gz |\
 cut -f1,2  >> $2/Filtred/SNP_2exclude.tsv
echo "non PASS SNP : "$(bcftools view -H -f.,DP_min_depth,QUAL_quality,FS_fisherstrand,QD_quality_by_depth,SOR_strand_odds_ratio,high_missing,high_heterozygosity $2/Filtred/combined.final.vcf.gz | wc -l) >> $2/Filtred/summary.txt


# Write singleton 
vcftools --gzvcf $2/Filtred/combined.final.vcf.gz --singletons --out $2/Filtred/to_remove

sed '1d' $2/Filtred/to_remove.singletons | awk '{print $1":"$2}' >> $2/Filtred/SNP_2exclude.txt
sed '1d' $2/Filtred/to_remove.singletons | cut -f1,2 >> $2/Filtred/SNP_2exclude.tsv
echo "doubletons : "$(sed '1d' $2/Filtred/to_remove.singletons | wc -l) >> $2/Filtred/summary.txt

## Hard filtering

bcftools view -T ^$2/Filtred/SNP_2exclude.tsv $2/Filtred/combined.final.vcf.gz -O z > $2/Filtred/combined.final.hard_filter.vcf.gz

# Store end time
end=$(date +%s)
# Calculate running time
runtime=$((end - start))
# Echo running time
echo "Script executed successfully. Running time: $runtime seconds."
