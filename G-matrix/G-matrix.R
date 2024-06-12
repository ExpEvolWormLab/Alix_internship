## Script to compute G-matrix using MCMCglmm 
# Same as FM script

# To change
args <- commandArgs(trailingOnly = TRUE)
M_file <- args[1] 
condition <- args[2]
output <- args[3]

M_file <- 'A6_t_convert_genotype.csv'
setwd("~/Documents/Worms/GBLUP")
condition <- 'NGM'

# Load packages
options(rgl.useNULL=TRUE)
library(data.table)
library(lattice)
library(MCMCglmm)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)

M <- as.data.frame(fread(M_file))
rownames(M) <- M[,1]
M[,1] <- NULL
M <- as.matrix(M)


# Load pheno data
pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
pheno_subset <- pheno[,c('pop_label','temperature','rel_humidity',"session_id",
                  'logD','date_str','env_label',
                  "SF",
                  "SB",
                  "FS",
                  "FB",
                  "BS",
                  "BF")]

## Kept only lines with phenotypes
pheno_subset<-pheno_subset[pheno_subset$env_label==condition,]
M <- M[rownames(M) %in% pheno_subset$pop_label,]


# Scale covariates and reduce traits
vect_P_traits <- c("SF",
                   "SB",
                   "FS",
                   "FB",
                   "BS",
                   "BF")

pheno_subset$temperature <- as.numeric(pheno_subset$temperature)
pheno_subset$rel_humidity <- as.numeric(pheno_subset$rel_humidity)
pheno_subset$logD <- as.numeric(pheno_subset$logD)
for(j in c('temperature',"rel_humidity","logD")){
  pheno_subset[,j][is.na(pheno_subset[,j])] <- mean(pheno_subset[,j], na.rm = TRUE)
  pheno_subset[,j] <- (pheno_subset[,j]-mean(pheno_subset[,j], na.rm = TRUE))/sd(pheno_subset[,j], na.rm = TRUE)
}
for(j in vect_P_traits){
  pheno_subset[,j]<-(pheno_subset[,j]-mean(pheno_subset[,j], na.rm = TRUE))
}



#### G-matrix

# Define prior model for GBLUP
nb_trait = length(vect_P_traits)
phen.var = diag(nb_trait) * diag(var(subset(pheno_subset, select = vect_P_traits)))
prior_mod <- list(G = list(G1 = list(V = phen.var/3, nu = nb_trait),
                           G2 = list(V = phen.var/3, nu = nb_trait)), 
                  R = list(V = phen.var/3, nu = nb_trait))
# MCMCglmm

if (length(unique(pheno_subset$session_id)) != 1) {
  model_MCMC_WI <- MCMCglmm(cbind(c(SF, SB, FS, FB, BS, BF)) ~ (logD + rel_humidity + temperature + session_id)^4 + trait - 1, 
                            random = ~us(trait):pop_label + us(trait):date_str, 
                            rcov = ~us(trait):units, 
                            family = rep("gaussian", nb_trait), 
                            prior = prior_mod, 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100)
} else {
  model_MCMC_WI <- MCMCglmm(cbind(c(SF, SB, FS, FB, BS, BF)) ~ (logD + rel_humidity + temperature)^4 + trait - 1, 
                            random = ~us(trait):pop_label + us(trait):date_str, 
                            rcov = ~us(trait):units, 
                            family = rep("gaussian", nb_trait), 
                            prior = prior_mod, 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100)
}



# # Save model results
write.csv(model_MCMC_WI$Sol, gsub('XXX', output, 'XXX_MCMCmodel_Sol.csv'), quote = FALSE, row.names = FALSE)
write.csv(model_MCMC_WI$VCV, gsub('XXX', output, 'XXX_MCMCmodel_VCV.csv'), quote = FALSE, row.names = FALSE)
write.csv( posterior.mode(model_MCMC_WI$VCV), gsub('XXX', output, 'XXX_MCMCmodel_VCV_posterior.mode.csv'), quote = FALSE, row.names = FALSE)

# Print model summary
print(summary(model_MCMC_WI))

post_dist = posterior.mode(model_MCMC_WI$VCV)
WI_G_mat=matrix(post_dist[1:nb_trait^2],nb_trait, nb_trait)
write.csv( WI_G_mat, gsub('XXX', output, 'XXX_MCMCmodel_VCV_G_matrix.csv'), quote = FALSE, row.names = FALSE)


# # Exit with a status code
# if (length(unique(pheno_subset$session_id)) != 1) {
#   quit(status = 0)
# } else {
#   quit(status = 1)
# }
