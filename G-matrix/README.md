# Set of Scripts Dedicated to G-Matrix Computation and Analysis

## G-matrix.R
### Terminal Script
```bash
Rscript G-matrix.R M_file pheno_file condition output
```
Use MCMCglmm to fit a BLUP model as done previously in FM work.

### Inputs:
- **M_file**: Path to the transposed genotype matrix.
- **pheno_file**: Path to the phenotypes file.
- **condition**: Environment condition (e.g., NGM/NaCl).
- **output**: Base name for the output files.

### Outputs:
- **MCMCmodel_Sol.csv**: Distribution a posteriori of breeding value for each effect for each line.
- **MCMCmodel_VCV.csv**: Variance-Covariance a posteriori for each effect for each line.
- **MCMCmodel_VCV_G_matrix.csv**: Matrix of variance covariance between our six phenotypes.

---

## Get_G_matrix.R
### Local Script
Get the G matrix from the MCMCmodel_VCV.csv computed during multivariate and univariate analysis.

### Input:
- **MCMCmodel_VCV.csv**

### Output:
- **G_matrix.csv**

---

## Compare_G.R
### Local Script
Compare different G_matrix.csv files. Just need to change name to seek files corresponding to name_G_matrix.csv.

### Output:
- Plot of correlation between G_matrix.
