#!/bin/bash

## Script to summarize number of SNPs by population afetr imputation and pruning
## Need to be launch in the same directory than launch_imputation_pruning.sh 

echo -e 'Pop\tSNPs_beforeI\tnSNPs_afterI\tnSNPS_afterF\tP_0.99\tP_0.9\tP_0.8\tP_0.7\tP_0.6\tP_0.5\tP_0.4\tP_0.3\tP_0.2\tP_0.1' > summary_imputation.txt

for i in A6 CA EEV GA GM GT LR SMR
do
	pop=${i%%_*}
	echo $pop
	nSNPs_b=$(bcftools view -H "$pop"_final.hard.vcf.gz | wc -l)
    	nSNPs=$(bcftools view -H imputed_All_"$pop".vcf.gz | wc -l)
	nSNPS_afterF=$(bcftools view -H "$i"/Filtred/"$pop".imputed.SNP.filtred.final.vcf.gz | wc -l)
    	l="$pop\t$nSNPs_b\t$nSNPs\t$nSNPS_afterF\t"
    	for k in 0.99 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1
	do
        	n=$(bcftools view -H "$i"/Filtred/Pruned/pruned."$k".vcf.gz | wc -l)
        	l+="\t$n"
    	done
    	echo -e "$l" >> summary_imputation.txt
done
