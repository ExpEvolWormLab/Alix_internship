## Script to develop GWAS using MCMCglmm: GBLUP step

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
matrix <- args[1]  # Path to the matrix file
M_file <- args[2]  # Path to the M file
condition <- args[3]  # Experimental condition
output <- args[4]  # Output directory

# Load necessary packages
options(rgl.useNULL=TRUE)
library(data.table)
library(lattice)
library(MCMCglmm)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)

# Load and process the kinship matrix
Ginv <- as.matrix(read.csv(matrix))
colnames(Ginv) <- gsub('CeMee', '', colnames(Ginv))  # Clean column names
colnames(Ginv) <- gsub('_sorted', '', colnames(Ginv))  # Clean column names
rownames(Ginv) <- colnames(Ginv)  # Set row names to match column names

# Load and process the M matrix
M <- as.data.frame(fread(M_file))
rownames(M) <- M[, 1]  # Set row names
M[, 1] <- NULL  # Remove the first column
M <- as.matrix(M)  # Convert to matrix

print(paste0(output, ': matrices loaded'))

# Load phenotype data
pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
pheno <- pheno[, c('pop_label', 'temperature', 'rel_humidity', "session_id", 'logD', 'date_str', 'env_label', "SF", "SB", "FS", "FB", "BS", "BF")]

# Keep only lines with phenotypes
pheno_subset <- pheno[pheno$pop_label %in% colnames(Ginv), ]
pheno_subset <- pheno_subset[pheno_subset$env_label == condition, ]
Ginv <- Ginv[rownames(Ginv) %in% pheno_subset$pop_label, rownames(Ginv) %in% pheno_subset$pop_label]
M <- M[rownames(M) %in% pheno_subset$pop_label, ]

# Scale covariates and standardize traits
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")
pheno_subset$temperature <- as.numeric(pheno_subset$temperature)
pheno_subset$rel_humidity <- as.numeric(pheno_subset$rel_humidity)
pheno_subset$logD <- as.numeric(pheno_subset$logD)

# Standardize covariates
for (j in c('temperature', "rel_humidity", "logD")) {
  pheno_subset[, j][is.na(pheno_subset[, j])] <- mean(pheno_subset[, j], na.rm = TRUE)
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE)) / sd(pheno_subset[, j], na.rm = TRUE)
}

# Center traits
for (j in vect_P_traits) {
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE))
}

# Plot phenotype data
plots <- lapply(c(vect_P_traits, 'temperature', 'rel_humidity', 'logD'), function(i) {
  ggplot(pheno_subset) + 
    geom_histogram(aes_string(x = i), fill = "grey", color = "black", bins = floor(nrow(pheno_subset) / 3)) +
    theme_bw()
})
pdf(gsub('XXX', output, "pheno_subset_plots_XXX.pdf"))
grid.arrange(grobs = plots, ncol = 3)
dev.off()

#### GBLUP ####
# Model as: y = Xb + Wr + u + e 
# X: fixed effects
# W: random effects
# u: follows normal distribution with variance of A*Var(u)

# Define prior model for GBLUP
nb_trait = length(vect_P_traits)
phen.var = diag(nb_trait) * diag(var(subset(pheno_subset, select = vect_P_traits)))
prior_mod <- list(G = list(G1 = list(V = phen.var / 3, nu = nb_trait),
                           G2 = list(V = phen.var / 3, nu = nb_trait)), 
                  R = list(V = phen.var / 3, nu = nb_trait))

# Convert Ginv to sparse matrix format
Ginv <- as(Ginv, "dgCMatrix")

# Perform GBLUP using MCMCglmm
if (length(unique(pheno_subset$session_id)) != 1) {
  model_MCMC_WI <- MCMCglmm(cbind(SF, SB, FS, FB, BS, BF) ~ (logD + rel_humidity + temperature + session_id) + trait - 1, 
                            random = ~us(trait):pop_label + us(trait):date_str, 
                            rcov = ~us(trait):units, 
                            family = rep("gaussian", nb_trait), 
                            prior = prior_mod, 
                            ginverse = list(pop_label = Ginv), 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100, 
                            pr = TRUE)
} else {
  model_MCMC_WI <- MCMCglmm(cbind(SF, SB, FS, FB, BS, BF) ~ (logD + rel_humidity + temperature) + trait - 1, 
                            random = ~us(trait):pop_label + us(trait):date_str, 
                            rcov = ~us(trait):units, 
                            family = rep("gaussian", nb_trait), 
                            prior = prior_mod, 
                            ginverse = list(pop_label = Ginv), 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100, 
                            pr = TRUE)
}

# Save model results
write.csv(model_MCMC_WI$Sol, gsub('XXX', output, 'XXX_MCMCmodel_Sol.csv'), quote = FALSE, row.names = FALSE)
write.csv(model_MCMC_WI$VCV, gsub('XXX', output, 'XXX_MCMCmodel_VCV.csv'), quote = FALSE, row.names = FALSE)

# Print model summary
print(summary(model_MCMC_WI))
