## Script to backsolve multivariate analysis ##

args <- commandArgs(trailingOnly = TRUE)
file <- args[1] # Sol
file1 <- args[2] # M
file2 <- args[3] # snp_pos

# Name of output
output <- args[4]

library(data.table)

model_MCMC_WI_Sol <- read.csv(file)
M <- as.data.frame(fread(file1))
rownames(M) <- M[,1]
M[,1] <- NULL
M <- as.matrix(M)
snp_positions <- read.csv(file2)

vect_P_traits <- c("SF",
                   "SB",
                   "FS",
                   "FB",
                   "BS",
                   "BF")

for(p in vect_P_traits){
  # get the vector of breeding value g for each trait 
  col_subset <- colnames(model_MCMC_WI_Sol)[grepl(paste0(p,'.pop_label'),colnames(model_MCMC_WI_Sol))]
  breeding_values <- data.frame(model_MCMC_WI_Sol[,col_subset])
  new_colnames <- gsub(".*pop_label\\.", "", colnames(breeding_values))
  colnames(breeding_values)<-new_colnames
  M <- M[rownames(M) %in% new_colnames,]
  
  # Backsolving
  # Calculer la transposée de M
  M_transpose <- t(M)
  
  # Calculer le produit MM'
  MM_prime <- M %*% M_transpose
  
  # Calculer l'inverse de MM'
  MM_prime_inverse <- solve(MM_prime)
  
  breeding_values_matrice <- as.matrix(breeding_values)
  
  SNPs_effects <-  M_transpose %*% MM_prime_inverse %*% t(breeding_values_matrice)
  
  write.csv(data.frame(cbind(snp_positions,SNPs_effects)),file = gsub('LLL',output,gsub('XXX',p,'SNPs_effects_XXX_LLL.csv')),quote = FALSE)
  write.csv(breeding_values,file = gsub('LLL',output,gsub('XXX',p,'BreedingValues_XXX_LLL.csv')),quote = FALSE)
}
