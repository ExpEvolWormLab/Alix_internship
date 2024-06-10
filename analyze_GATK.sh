#!/bin/bash
# First argument : file to analyse
# Second name of output

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
