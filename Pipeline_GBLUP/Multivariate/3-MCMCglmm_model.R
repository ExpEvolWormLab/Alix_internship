## Script to developp GWAS using MCMCglmm : GBLUP step


args <- commandArgs(trailingOnly = TRUE)
matrix <- args[1]
M_file <- args[2] 
condition <- args[3]
output <- args[4]

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

print(paste0(output,': matrices loaded'))
# Load pheno data
pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
pheno <- pheno[,c('pop_label','temperature','rel_humidity',"session_id",
                  'logD','date_str','env_label',
                  "SF",
                  "SB",
                  "FS",
                  "FB",
                  "BS",
                  "BF")]

## Kept only lines with phenotypes
pheno_subset<-pheno[pheno$pop_label %in% colnames(Ginv),]
pheno_subset<-pheno_subset[pheno_subset$env_label==condition,]
#pheno_subset<-pheno_subset[complete.cases(pheno_subset[,c('temperature', 'rel_humidity', 'logD')]),]
Ginv <- Ginv[rownames(Ginv) %in% pheno_subset$pop_label,rownames(Ginv) %in% pheno_subset$pop_label]
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


# Plot phenotype data
plots <- lapply(c(vect_P_traits, 'temperature', 'rel_humidity', 'logD'), function(i) {
  ggplot(pheno_subset) + 
    geom_histogram(aes_string(x = i), fill = "grey", color = "black", bins = floor(nrow(pheno_subset) / 3)) +
    theme_bw()
})
pdf(gsub('XXX',output,"pheno_subset_plots_XXX.pdf"))
grid.arrange(grobs = plots, ncol = 3)
dev.off()


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

# # Save model results
write.csv(model_MCMC_WI$Sol, gsub('XXX', output, 'XXX_MCMCmodel_Sol.csv'), quote = FALSE, row.names = FALSE)
write.csv(model_MCMC_WI$VCV, gsub('XXX', output, 'XXX_MCMCmodel_VCV.csv'), quote = FALSE, row.names = FALSE)

# Print model summary
print(summary(model_MCMC_WI))

# # Exit with a status code
# if (length(unique(pheno_subset$session_id)) != 1) {
#   quit(status = 0)
# } else {
#   quit(status = 1)
# }
