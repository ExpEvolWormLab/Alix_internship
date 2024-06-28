# Automated GBLUP Analysis Procedure

## Running the Analysis

### Multivariate Analysis:
```bash
bash launch_GBLUP_pop.sh file.tsv
```

### Univariate Analysis:
```bash
bash launch_GBLUP_pop.Univariate.sh file.tsv
```

**`file.tsv`**: A TSV file where the first column contains the path of the VCF file to use for a population, and the second column contains the name of the population.

## Details

### `launch_GBLUP_pop.sh`

- Launches `launch_GBLUP.sh`

### `launch_GBLUP.sh`

**`1-Kinship.R`:**
   - Computes the kinship matrix in three different ways (Hoffman, NOIA, VanRaden)
   - Outputs:
     - `*_snp_positions.csv`: SNP positions
     - `*_t_convert_genotype.csv`: Transposed genotype matrix
     - `Kinship_matrix_VanRaden.csv`: VanRaden kinship matrix
     - `Kinship_matrix_Luke.csv`: Hoffman kinship matrix
     - `Kinship_matrix_noia.csv`: Noia kinship matrix

**Launches `launch_GBLUP.parallele.sh`:**
   - Proceeds with the rest of the procedure using the three ways of computing GRM

### `launch_GBLUP.parallele.sh`

**`2-Invert_GRM.R`:**
   - Computes the inverse of the kinship matrix if it's ill-conditioned or not positive definite, and modifies it slightly
   - Output:
     - `Inverted_kinship_matrix.csv`: Inverted kinship matrix

**`3-MCMCglmm_model.R`:**
   - Fits the GBLUP model: \( Y = Xb + Wr + Zu + e \)
     - \( Y \): Phenotypes
     - \( X \): Matrix of fixed effects (temperature, density, humidity, session (only if not always the same))
     - \( W \): Random effects (block effect)
     - \( Z \): Line identity (variance constrained by the kinship matrix)
     - \( e \): Error term following a normal distribution centered at 0
   - Outputs:
     - `MCMCmodel_Sol.csv`: A posteriori distribution for each line for each trait and effect
     - `MCMCmodel_VCV.csv`: Variance-covariance between traits and effects

**`4-Diagnostic.R`:**
   - Provides a visual diagnostic of the convergence of the MCMCglmm model
   - Output:
     - `Model_pdf_MCMC_autocorrelation.pdf`: Plots showing the autocorrelation per iteration

**`5-Backsolving.R`:**
   - Uses `MCMCmodel_Sol.csv` to get SNP effects
   - Outputs:
     - `BreedingValues.csv`: Table storing breeding value distribution a posteriori for each line for each trait
     - `SNPs_effects.csv`: Table storing SNP effects a posteriori distribution for each line for each trait
    
**`5-Backsolving_Pruned_Data.R`:**
   -  Uses `MCMCmodel_Sol.csv` and M matrix from pruned data to get SNP effects
   -  Outputs:
     - `BreedingValues.csv`: Table storing breeding value distribution a posteriori for each line for each trait
     - `SNPs_effects.csv`: Table storing SNP effects a posteriori distribution for each line for each trait
