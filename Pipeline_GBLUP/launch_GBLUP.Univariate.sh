#!bin/bash

## Script which launch each steps of the GBLUP and Backsolving pipeline
## First argument name of vcf to use
## Second argument name of the output (will be add on the extension of created files, all files will be created where this script is launch)

## Rscripts have to be in the same directory

Rscript 1-Kinship.Univariate.R $1 $2
echo '$2: Kinship matrix executed'


for i in "VanRaden_" "noia_" "Luke_"
do
	bash launch_GBLUP.parallele.Univariate.sh "Kinship_matrix_${i}$2.csv" "${i}$2" "$2_t_convert_genotype.csv" "$2_snp_positions.csv"
done
wait

mkdir -p $2

mv *$2* $2
