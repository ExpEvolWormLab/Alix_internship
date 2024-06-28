## Script to backsolve with pruned matrix for UNIVARIATE##
output <- 'Uni_VanRaden_A6_NaCl_0.99'
file1 <- 'pruned.0.99.vcf.gz'
condition <- 'NaCl'


vect_P_traits <- c("SF",
                   "SB",
                   "FS",
                   "FB",
                   "BS",
                   "BF")

library(data.table)
library(vcfR)


setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")

# Read vcf file
vcf <- read.vcfR(file1)

# Extract genotype information
genotype_info <- vcf@gt

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]

# remove FORMAT column if it exists
if ("FORMAT" %in% colnames(genotype_info)) {
  genotype_info <- subset(genotype_info, select = -c(FORMAT))
}

# Function to convert genotypes
get_genotype <- function(genotype) {
  ifelse(substr(genotype, 1, 3) %in% c("0/0", "0|0"), 0,
         ifelse(substr(genotype, 1, 3) %in% c("1/1", "1|1"), 2,
                ifelse(substr(genotype, 1, 3) %in% c("./.", ".|."), NaN,NaN)))
}

convert_genotype <- apply(genotype_info, 2, get_genotype)

snp_positions <- as.data.frame(snp_positions[complete.cases(convert_genotype),])
convert_genotype <- convert_genotype[complete.cases(convert_genotype),]

# Function to clean the column names
clean_colnames <- function(colnames) {
  sapply(colnames, function(name) {
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) > 1 && parts[1] == parts[2]) {
      return(parts[1])
    } else {
      return(name)
    }
  })
}

# Clean the column and row names
colnames(convert_genotype) <- clean_colnames(colnames(convert_genotype))
M <- as.matrix(t(convert_genotype))
rownames(M) <- gsub('CeMee','',rownames(M))
rownames(M) <- gsub('_sorted','',rownames(M))


for(p in vect_P_traits){
  setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
  file <- gsub('condition',condition,gsub('TTT', p,'VanRaden_A6_condition_TTT_MCMCmodel_Sol.csv'))
  model_MCMC_WI_Sol <- read.csv(file)


  # get the vector of breeding value g for each trait 
  col_subset <- colnames(model_MCMC_WI_Sol)[grepl('pop_label',colnames(model_MCMC_WI_Sol))]
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
  
  setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")
  write.csv(data.frame(cbind(snp_positions,SNPs_effects)),file = gsub('LLL',output,gsub('XXX',p,'SNPs_effects_XXX_LLL.csv')),quote = FALSE)
  write.csv(breeding_values,file = gsub('LLL',output,gsub('XXX',p,'BreedingValues_XXX_LLL.csv')),quote = FALSE)
}
