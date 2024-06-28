## Script to compute G-matrix using MCMCglmm

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
M_file <- args[1]  # Path to the M file
pheno_file <- args[2]
condition <- args[3]  # Condition to filter phenotype data
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

# Load and process the M matrix
M <- as.data.frame(fread(M_file))
rownames(M) <- M[, 1]  # Set row names
M[, 1] <- NULL  # Remove the first column
M <- as.matrix(M)  # Convert to matrix

# Load phenotype data
pheno <- read.csv(pheno_file)
pheno_subset <- pheno[, c('pop_label', 'temperature', 'rel_humidity', "session_id",
                          'logD', 'date_str', 'env_label', "SF", "SB", "FS", "FB", "BS", "BF")]

# Keep only lines with phenotypes for the specified condition
pheno_subset <- pheno_subset[pheno_subset$env_label == condition, ]
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

#### G-matrix

# Define prior model for GBLUP
nb_trait = length(vect_P_traits)
phen.var = diag(nb_trait) * diag(var(subset(pheno_subset, select = vect_P_traits)))
prior_mod <- list(G = list(G1 = list(V = phen.var / 3, nu = nb_trait),
                           G2 = list(V = phen.var / 3, nu = nb_trait)), 
                  R = list(V = phen.var / 3, nu = nb_trait))

# Fit the MCMCglmm model
if (length(unique(pheno_subset$session_id)) != 1) {
  model_MCMC_WI <- MCMCglmm(cbind(SF, SB, FS, FB, BS, BF) ~ (logD + rel_humidity + temperature + session_id)^4 + trait - 1, 
                            random = ~us(trait):pop_label + us(trait):date_str, 
                            rcov = ~us(trait):units, 
                            family = rep("gaussian", nb_trait), 
                            prior = prior_mod, 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100)
} else {
  model_MCMC_WI <- MCMCglmm(cbind(SF, SB, FS, FB, BS, BF) ~ (logD + rel_humidity + temperature)^4 + trait - 1, 
                            random = ~us(trait):pop_label + us(trait):date_str, 
                            rcov = ~us(trait):units, 
                            family = rep("gaussian", nb_trait), 
                            prior = prior_mod, 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100)
}

# Save model results
write.csv(model_MCMC_WI$Sol, gsub('XXX', output, 'XXX_MCMCmodel_Sol.csv'), quote = FALSE, row.names = FALSE)
write.csv(model_MCMC_WI$VCV, gsub('XXX', output, 'XXX_MCMCmodel_VCV.csv'), quote = FALSE, row.names = FALSE)

# Print model summary
print(summary(model_MCMC_WI))

# Compute and save the G-matrix
post_dist = posterior.mode(model_MCMC_WI$VCV)
WI_G_mat = matrix(post_dist[1:nb_trait^2], nb_trait, nb_trait)
write.csv(WI_G_mat, gsub('XXX', output, 'XXX_MCMCmodel_VCV_G_matrix.csv'), quote = FALSE, row.names = FALSE)
