#!/bin/bash

#Script to do pruning with different r2
# feature 1 : file to prune (.vcf format)
# feature 2 : where to store results


#Create directory where the results will be stored
mkdir -p $2/Pruned

# create chrom map 
#bcftools view -H $1 | cut -f 1 | uniq | awk '{print $0"\t"$0}' > $2/Pruned/chrom-map.txt

# convert vcf to plink format
#vcftools --gzvcf $1 --plink --out $2/Pruned/plink.data --chrom-map $2/Pruned/chrom-map.txt

#bcftools annotate --set-id +'%CHROM\_%POS' $i -O z > 

# pruning
for r2 in 0.99 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
do
	plink --vcf $1 --indep-pairwise 200 10 $r2 --out $2/Pruned/pruned.$r2 --allow-extra-chr --double-id
	plink --vcf $1 --extract $2/Pruned/pruned.$r2.prune.in --recode vcf --out $2/Pruned/pruned.$r2 --allow-extra-chr --double-id # --update-name --set-missing-var-ids @:# 
	bgzip $2/Pruned/pruned.$r2.vcf
done
