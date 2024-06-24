#!/bin/bash

# first feature : final file of imputation
# second : name of the directory where the results'll be store

# To change :
path_gatk=$3

set -e #Stop if error
start=$(date +%s) #store time

# Create directories
mkdir -p $2/Filtred

# Variant Filtration
# Keep only SNPs

## Sometimes SelectVariant keep ref allele with only * : doesn't account as a SNP, remove them with :
bcftools view -v snps $1 -O z > $2/Filtred/imputed.SNP.vcf.gz
$path_gatk IndexFeatureFile -I $2/Filtred/imputed.SNP.vcf.gz

# Spot heterozygotes
$path_gatk VariantFiltration \
	-V $2/Filtred/imputed.SNP.vcf.gz \
	--genotype-filter-expression "isHet == 1.0"	--genotype-filter-name "is_het" \
	-O $2/Filtred/imputed.SNP.het.vcf.gz

# Apply high heterozygosity filter
bcftools filter --soft-filter='high_heterozygosity' --mode + --include 'F_PASS(GT="het")<=0.15' -O z  $2/Filtred/imputed.SNP.het.vcf.gz > $2/Filtred/imputed.SNP.het1.vcf.gz
$path_gatk IndexFeatureFile -I $2/Filtred/imputed.SNP.het1.vcf.gz

#Put ./. ie missing value for  heterozygotes genotype
$path_gatk SelectVariants \
	-V $2/Filtred/imputed.SNP.het1.vcf.gz \
	--set-filtered-gt-to-nocall \
	-O $2/Filtred/imputed.SNP.without_hetero.vcf.gz

# Apply high missing filter
bcftools filter --soft-filter='high_missing' --mode + --include 'F_MISSING<=0.95' $2/Filtred/imputed.SNP.without_hetero.vcf.gz -O z >  $2/Filtred/imputed.SNP.filtred.vcf.gz


# Hard filtering, remove sites : nonPass, multiallelic, monomorphic
bcftools view -m2 -M2 -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' -f PASS $2/Filtred/imputed.SNP.filtred.vcf.gz -O z > $2/Filtred/imputed.SNP.hard_filtred.vcf.gz


vcftools --maf 0.05 --gzvcf $2/Filtred/imputed.SNP.hard_filtred.vcf.gz --out $2/Filtred/imputed.SNP.filtred.final.vcf --recode
mv $2/Filtred/imputed.SNP.filtred.final.vcf.recode.vcf $2/Filtred/imputed.SNP.filtred.final.vcf
bgzip $2/Filtred/imputed.SNP.filtred.final.vcf

# Store end time
end=$(date +%s)
# Calculate running time
runtime=$((end - start))
# Echo running time
echo "Script executed successfully. Running time: $runtime seconds."
