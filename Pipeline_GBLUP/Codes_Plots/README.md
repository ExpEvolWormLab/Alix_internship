# Set of Scripts to Plot GWAS by GBLUP Results

## Breeding Values

### Univariate

#### Compare_BV.Pruning.Univariate.R
Script to compare breeding values of pruned and non-pruned data.

**Inputs:** 
- `name`: Pattern to search for pruned MCMCmodel_Sol.csv files.
- `name1`: Pattern to search for non-pruned MCMCmodel_Sol.csv files.
- `M_file`: Converted genotype file.
- `pheno_file`: File containing phenotypic measures.

**Output:** 
- `Median_BV_among_pruning.pdf`: Plot of median and standard deviation of breeding values among pruned and non-pruned steps.

#### Plot_BV_Univariate.R
Script to plot results of breeding values for univariate analysis.

**Inputs:** 
- `output`: Output name.
- `M_file`: Converted genotype file.
- `pheno_file`: File containing phenotypic measures.

**Outputs:** 
- `Trait_distribution_each_lines*pdf`: Breeding value distribution (median and SD) for each line along the trait.
- `forest_plot*pdf`: Depicts breeding values that are credible (95% credibility interval doesn't overlap 0).
- `Trait_VS_traits_plot*pdf`: Plots breeding values for each trait against each other, computes Pearson coefficients to attest the correlation.

### Multivariate

#### Plot_BV.R
Script to plot results of breeding values for multivariate analysis.

**Inputs:** 
- `output`: Output name.
- `M_file`: Converted genotype file.
- `pheno_file`: File containing phenotypic measures.

**Outputs:** 
- `Trait_distribution_each_lines*pdf`: Breeding value distribution (median and SD) for each line along the trait.
- `forest_plot*pdf`: Depicts breeding values that are credible (95% credibility interval doesn't overlap 0).
- `Trait_VS_traits_plot*pdf`: Plots breeding values for each trait against each other, computes Pearson coefficients to attest the correlation.

## SNP Effects

### First Step
Running:

```bash
Rscript Summarize_SNPsEffect.R SNPs_effect output_name
```

Or, if considering pruned data:

```bash
Rscript Summarize_SNPsEffect_pruned.R SNPs_effect output_name
```

**Output:** 
- `Summary*csv`: Name of SNP (CHROM_POS), CHROM, Median (of SNPs effect a posteriori distribution), lower (of credibility interval), upper, Credibility.

### Univariate

#### Compare_Credible_SNPs_Univariate.R
**Inputs:** 
- `output`: Pattern to search for summary files and to name output files.
- `output_dir`: Path where results are going to be stored.

**Outputs:** 
- `Shared_credible_SNPs*tsv`: File with two columns (Comparison, Nbr_SNPs), comparing each trait and storing the number of credible SNPs they share.
- `Shared_credible_SNPs*pdf`: Plot of the previous files.
- `Repartition_SNPs_along_K`: For each SNP, plot it along the chromosome in function of the number of groups they are associated with.

#### SNPs_effects_Univariate.R
Script to analyze SNP effects for univariate analysis.

**Inputs:** 
- `output`: Pattern to search for summary files and to name output files.
- `SNP_Count`: File which stores for SNPs the number of groups they are associated with.

**Outputs:** 
- `Manhattan Plot`: Median of posterior SNP effects along the chromosome colored by credibility.
- `Density Plot`: Density function of the median of the posterior distribution of SNP effects.

### Multivariate

#### Compare_Credible_SNPs.R
**Inputs:** 
- `output`: Pattern to search for summary files and to name output files.
- `output_dir`: Path where results are going to be stored.

**Outputs:** 
- `Shared_credible_SNPs*tsv`: File with two columns (Comparison, Nbr_SNPs), comparing each trait and storing the number of credible SNPs they share.
- `Shared_credible_SNPs*pdf`: Plot of the previous files.
- `SNP_Count`: File which stores for SNPs the number of groups they are associated with.
- `Repartition_SNPs_along_K`: For each SNP, plot it along the chromosome in function of the number of groups they are associated with.

#### SNPs_effects.R
Script to analyze SNP effects for multivariate analysis.

**Inputs:** 
- `output`: Pattern to search for summary files and to name output files.
- `SNP_Count`: File which stores for SNPs the number of groups they are associated with.

**Outputs:** 
- `Manhattan Plot`: Median of posterior SNP effects along the chromosome colored by credibility.
- `Density Plot`: Density function of the median of the posterior distribution of SNP effects.

#### MP_across_conditions.R
Script to compare credible SNPs between conditions (NaCl/NGM).

**Inputs:** 
- `output`: Pattern to search for summary files and name of the output.
- `pruning`: Pruning pattern to search for summary files.

**Output:** 
- `Manhattan Plot`: Median posterior effect of SNPs along the chromosome colored by association in the two environments.
