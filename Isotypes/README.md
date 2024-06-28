# Set of Scripts to Obtain Isotype Groups

This document provides an overview of the scripts used to obtain isotype groups. The steps include running the NaN distribution analysis, generating and correcting concordance matrices, identifying problematic isotypes, and converting the phenotype table. Each script's usage and output files are described to facilitate the process of obtaining and analyzing isotype groups.

## 1. NaN Distribution

### Script: `1-NaN_Distribution.R`

#### Usage:
- This script should use a hard-filtered VCF file to consider only PASS SNPs.

#### Outputs:
- `Na_distribution.tsv`: Tab-separated values file containing the NaN distribution data.
- `histo_NA.pdf`: Histogram of the NaN distribution used to set a cutoff.

## 2. Concordance Matrix

### Script: `2-Concordance_matrix.R`

#### Usage:
- This script should use a hard-filtered VCF file to consider only PASS SNPs.
- Concordance is calculated as the number of SNPs with the same allele in both lines divided by the number of SNPs observed in both lines.

#### Outputs:
- `Doublons.csv`: Contains names of lines sequenced twice.
- `matrix_concordance.by_hand*.tsv`: Concordance matrix files.

## 3. Set Cutoffs for Isotypes

### Script: `3-Set_cutoffs_Isotypes.R`

#### Outputs:
- `Removed_line_concordance.csv`: Stores names of lines belonging to known isotypes with low concordance.
- `Removed_line_doublon.csv`: Stores names of lines considered twice by mistake.
- `matrix_concordance.by_hand.corrected.tsv`: Corrected concordance matrix, removing known isotypes inferior to the cutoff.

## 4. Identify Problems with Isotypes

### Script: `4-Pb_Isotypes.R`

#### Description:
- Computes isotype groups and identifies those shared between several groups.
- Computes a new concordance matrix.

#### Outputs:
- `matrix_concordance.by_hand.corrected1.tsv`: Corrected concordance matrix, excluding problematic lines.
- `Removed_line_isotype_groups.csv`: Stores names of lines shared between groups.

## 5. Isotypes

### Script: `5-Isotypes.R`

#### Outputs:
- `isotype_groups.csv`: Groups of isotypes.
- `Distribution_isotypes.csv`: Number of lines in each population belonging to an isotype group.
- `correspondance_isotypes.csv`: Lines belonging to an isotype group and the line representing it.
- `Removed_line_isotypes.csv`: Stores names of lines removed because they belong to an isotype group.


## Pipeline
![Results_GWAS_Worms](https://github.com/ExpEvolWormLab/Alix_internship/assets/83120878/552e0b51-3e97-4cc9-8eca-cddc3b6eb451)

## Convert Phenotype Table

### Script: `Convert_phenotype_table.R`

#### Description:
- Converts the phenotype table to modify the names of lines belonging to isotype groups as described in `correspondance_isotypes.csv`.

#### Output:
- `Final_Transition_rates_estimates_may2024_export.csv`: Converted phenotype table with updated line names.

