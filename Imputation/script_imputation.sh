#!/bin/bash

## Script which launch imputation by populations
# first feature name of the vcf file
# second feature name of the output directory

mkdir -p "$2"

for pop in A6 CA EEV GA GM GT LR SMR # For each population
do
    if [ ! -e "$2/${pop}_final.vcf.gz" ] # if the file by pop hasn't be compute yet
    then
        # Get the name of line belonging to the pop
        bcftools query -l "$1" | grep "$pop" > "$2/${pop}.sample.name"
        # Extract them from vcf file
        bcftools view "$1" -S "$2/${pop}.sample.name" -O z > "$2/${pop}_final.vcf.gz"
    fi
    wait
    for chrom in I II III
    do
        java -Xmx10g -jar beagle.22Jul22.46e.jar excludemarkers="SNP_2exclude.txt" gt="$2/${pop}_final.vcf.gz" out="$2/imputed_${chrom}_${pop}" chrom="$chrom" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 &
    done
    wait
    for chrom in IV V X
    do
        java -Xmx10g -jar beagle.22Jul22.46e.jar excludemarkers="SNP_2exclude.txt" gt="$2/${pop}_final.vcf.gz" out="$2/imputed_${chrom}_${pop}" chrom="$chrom" ne=1000 window=5 overlap=2 impute=true imp-segment=0.5 imp-step=0.01 cluster=0.0005 &
    done
    wait
done
