#!/bin/bash

## Script to analyse a vcf file 
## Need its output for Plot_stat_vcf.R 

# First argument : file to analyse
# Second name of output

# Outputs : pdf file with plots and summary, .vchk needed for Plot_stat_vcf.R 

bcftools query -l $1 > "sample.name"
bcftools stats -S "sample.name" $1 > "output_$2.vchk"
grep "^PSC" "output_$2.vchk" > "output_$2_PSC.vchk"
plot-vcfstats "output_$2.vchk" -p "output_$2"
cd "output_$2"
mv "summary.pdf" "summary_$2.pdf"
mv "summary_$2.pdf" ../.
cd ..
rm -r "output_$2"
rm "output_$2.vchk"
