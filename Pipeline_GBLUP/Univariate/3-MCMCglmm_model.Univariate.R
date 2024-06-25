## Script to developp GWAS using MCMCglmm : fit MCMCglmm model (GBLUP)

args <- commandArgs(trailingOnly = TRUE)
matrix <- args[1]
M_file <- args[2] 
condition <- args[3]
trait <- args[4]
output <- args[5]

matrix <- 'Inverted_kinship_matrix_VanRaden_GA.csv'
M_file <- 'GA_t_convert_genotype.csv'
condition <- 'NGM'
trait <- 'FS'

# Load packages
options(rgl.useNULL=TRUE)
library(data.table)
library(lattice)
library(MCMCglmm)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)

Ginv <- as.matrix(read.csv(matrix))
colnames(Ginv) <- gsub('CeMee','',colnames(Ginv))
colnames(Ginv) <- gsub('_sorted','',colnames(Ginv))
rownames(Ginv)<-colnames(Ginv)
M <- as.data.frame(fread(M_file))
rownames(M) <- M[,1]
M[,1] <- NULL
M <- as.matrix(M)

#print(paste0(output,': matrices loaded'))
# Load pheno data
pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv.txt')
pheno <- pheno[,c('pop_label','temperature','rel_humidity',"session_id",
                  'logD','date_str','env_label',
                  trait)]

## Kept only lines with phenotypes
pheno_subset<-pheno[pheno$pop_label %in% colnames(Ginv),]
pheno_subset<-pheno_subset[pheno_subset$env_label==condition,]
Ginv <- Ginv[rownames(Ginv) %in% pheno_subset$pop_label,rownames(Ginv) %in% pheno_subset$pop_label]
M <- M[rownames(M) %in% pheno_subset$pop_label,]


# Scale covariates and reduce traits
vect_P_traits <- trait

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

#### GBLUP ####
# Model as : y= Xb + Wr + u + e 
# X, fixed effects
# W, random effects
# u, follow normal distribution with variance of A*Var(u)

# Define prior model for GBLUP
nb_trait = length(vect_P_traits)
phen.var = diag(nb_trait) * diag(var(subset(pheno_subset, select = vect_P_traits)))
prior_mod <- list(G = list(G1 = list(V = phen.var/3, nu = nb_trait),
                           G2 = list(V = phen.var/3, nu = nb_trait)), 
                  R = list(V = phen.var/3, nu = nb_trait))


# Perform GBLUP using MCMCglmm
Ginv <- as(Ginv, "dgCMatrix")


if (length(unique(pheno_subset$session_id)) != 1) {
  fixed_effects <- paste(trait, "~ logD + rel_humidity + temperature + session_id - 1")
  formula = as.formula(paste(fixed_effects, collapse = " + "))
  model_MCMC_WI <- MCMCglmm(fixed=formula, 
                            random = ~pop_label + date_str, 
                            family = "gaussian", 
                            prior = prior_mod, 
                            ginverse = list(pop_label = Ginv), 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100, 
                            pr = TRUE)
} else {
  fixed_effects <- paste(trait, "~ logD + rel_humidity + temperature - 1")
  formula = as.formula(paste(fixed_effects, collapse = " + "))
  model_MCMC_WI <- MCMCglmm(fixed=formula, 
                            random = ~pop_label + date_str, 
                            family = "gaussian", 
                            prior = prior_mod, 
                            ginverse = list(pop_label = Ginv), 
                            data = pheno_subset, 
                            nitt = 110000, 
                            burnin = 10000, 
                            thin = 100, 
                            pr = TRUE)
}

write.csv(model_MCMC_WI$Sol,gsub('XXX',output,'XXX_MCMCmodel_Sol.csv'),quote = FALSE,row.names = FALSE)
write.csv(model_MCMC_WI$VCV,gsub('XXX',output,'XXX_MCMCmodel_VCV.csv'),quote = FALSE,row.names = FALSE)

summary.MCMCglmm(model_MCMC_WI)
