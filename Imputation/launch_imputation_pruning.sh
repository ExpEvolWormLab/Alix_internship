#!bin/bash
## Script which launch imputation by populations, filter it and prune it 
path_gatk='/mnt/data2/desmarais/gatk-4.5.0.0/gatk'


# First feature name of the vcf file
# Second feature name of the output directory

bash script_imputation.sh $1 $2

## Concat imputed file by chromosomes for each pop and launch filtration and pruning


for pop in A6 CA EEV GA GM GT LR SMR # For each population
do
	(
		bcftools concat --naive imputed_*_$pop.vcf.gz -O z > imputed_All_$pop.vcf.gz
  		# Name SNPs
		bcftools reheader -f chrom_lenght.fa.fai imputed_All_$pop.vcf.gz > imputed_All_$pop.H.vcf.gz
		bcftools annotate --set-id +'%CHROM\_%POS' imputed_All_$pop.H.vcf.gz  -O z > imputed_All_$pop.ID.vcf.gz
		bash script_filtration_imputation.sh imputed_All_$pop.ID.vcf.gz $pop
	) &
done

wait

bash launch_pruning.sh

wait
