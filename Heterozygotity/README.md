# Set of Scripts to Get Heterozygosity Along Chromosomes and to Plot It

This document provides an overview of how to get and plot heterozygosity along chromosomes using the provided scripts and tools. The steps include running the initial heterozygosity extraction script, followed by plotting the results locally using the provided R scripts to generate informative plots.

## 1. Get Heterozygosity Data

### Script: `1-Get_hetero.R`

#### Command to Run:
```bash
Rscript 1-Get_hetero.R het_flag.vcf
```

#### Parameters:
- `het_flag.vcf`: VCF file where heterozygote loci are flagged as `is_het` (e.g., `combined.annotated.filtered.vcf.gz`).
  - **Note**: Some parameters such as chromosomes to plot, populations to plot, or number of SNPs in a bin can be changed directly in the script.

#### Output:
- `Heterozygosity*.csv`: Contains the mean heterozygosity of lines for the SNP windows defined by the user along chromosomes.

## 2. Plot Heterozygosity Data

### Script: `2-Plotting_hetero.R`

#### Command to Run:
- This script should be run locally.

#### Requirements:
- Hyper divergent file is necessary.

#### Outputs:
- `Heterozygosity_HDregion.pdf`: Heterozygosity along chromosomes with hyper divergent regions plotted as black dots.
- `Heterozygosity_HDregion_Zoom_V.pdf`: Heterozygosity for the remarkable regions of chromosome V with hyper divergent regions plotted as black dots.
- `Randomisation.pdf`: Comparison of the observed proportion of heterozygotes associated with hyper divergent regions (HDR) for all chromosomes compared to the distribution obtained by chance.
- `Randomisation_Chrom5.pdf`: Comparison of the observed proportion of heterozygotes associated with HDR for the remarkable regions of chromosome V compared to the distribution obtained by chance.

## 3. Plot Heterozygosity for Chromosome 5

### Script: `Plotting_Chrom5.R`

#### Command to Run:
- This script should be run locally.

#### Requirements:
- Hyper divergent file is necessary.
- All population data should be in the working directory named as `Heterozygosity_pop.csv`.

#### Output:
- `Chromosome_5_Allpop.pdf`: Heterozygosity along chromosome V with hyper divergent regions plotted as black dots for all populations in the remarkable region of chromosome V.


