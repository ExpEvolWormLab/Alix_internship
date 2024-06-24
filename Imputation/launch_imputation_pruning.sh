#!bin/bash
## Script which launch imputation by populations, filter it and prune it 
path_gatk='/mnt/data2/desmarais/gatk-4.5.0.0/gatk'


# First feature name of the vcf file
# Second feature name of the output directory

bash script_imputation.sh $1 $2

## Concat imputed file by chromosomes for each pop and launch filtration and pruning


for pop in A6 CA EEV GA GM GT LR SMR # For each population
do
	bcftools concat imputed_*_$pop.vcf.gz -O z > imputed_All_$pop.vcf.gz
	bash script_filtration_imputation.sh imputed_All_$pop.vcf.gz $pop $path_gatk
done &

wait

bash launch_pruning.sh

wait
