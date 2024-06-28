# Scripts to Extract Population Data, Plot Environmental Variance, and Visualize Phenotypic Data

## Extract_pop.R
### Inputs
- **Final_lines.csv**: Line names that are retained after the line filtration step (944 lines).
- **Final_Transition_rates_estimates_may2024_export.csv**: File containing phenotypic measures for each line, corrected with isotype groups.
- **Output name**: Name of the output file.

### Outputs
- **output_name.csv**: Distribution of genotyped lines retained in the CeMee Panel.
- **output_name_phenotyped.csv**: Distribution of genotyped and phenotyped lines retained in the CeMee Panel.

## Plot_Environmental_Variance.R
### Inputs
- **Inverted_kinship_matrix_VanRaden_A6.csv**: Inverted kinship matrix.
- **Final_Transition_rates_estimates_may2024_export.csv**: File containing phenotypic measures for each line, corrected with isotype groups.
- **Output name**: Name of the output file.

### Outputs
- **Environmental_variance*pdf**: Plot depicting environmental variance for each line in both conditions (NaCl/NGM).

## Plot_subset_pheno.R
### Inputs
- **Inverted_kinship_matrix_VanRaden_A6.csv**: Inverted kinship matrix.
- **Final_Transition_rates_estimates_may2024_export.csv**: File containing phenotypic measures for each line, corrected with isotype groups.
- **Output name**: Name of the output file.

### Outputs
- **plots_pheno_subset.pdf**: Density function of phenotype values compared in both conditions for each trait.
- **Pheno_Trait_VS_traits_plot**: Plot depicting the correlation between traits using phenotypic measures and Pearson Coefficients.

## Plots_pheno_along_genotypes.R
### Description
Script to plot phenotype and genotype data as in Stephens and al. (2013).

### Inputs
- **Inverted_kinship_matrix_VanRaden_A6.csv**: Inverted kinship matrix.
- **Final_Transition_rates_estimates_may2024_export.csv**: File containing phenotypic measures for each line, corrected with isotype groups.
- **pruned.0.99.vcf.gz**: Pruned VCF file.
- **chrom**: Specific chromosome to plot.
- **pos**: Specific position to plot.
- **Output name**: Name of the output file.

### Outputs
- **Pheno_Geno_Trait_VS_traits_plot**: Plot depicting the correlation between FS and other traits using phenotypic measures and Pearson Coefficients, colored along the genotype of the specific position.
