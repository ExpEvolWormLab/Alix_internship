## Script to backsolve ##

# Load necessary libraries
library(data.table)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]  # Solution file (Sol)
file1 <- args[2]  # Genotype file (M)
file2 <- args[3]  # SNP positions file (snp_pos)
trait <- args[4]  # Trait to analyze

# Name of output file
output <- args[5]

# Load the solution data
model_MCMC_WI_Sol <- read.csv(file)

# Load the genotype data
M <- as.data.frame(fread(file1))
rownames(M) <- M[,1]  # Set row names to the first column
M[,1] <- NULL  # Remove the first column
M <- as.matrix(M)  # Convert to matrix

# Load the SNP positions data
snp_positions <- read.csv(file2)

# Define the trait vector
vect_P_traits <- trait

# Loop through each trait in the trait vector
for(p in vect_P_traits) {
  # Extract the vector of breeding values for the trait
  col_subset <- colnames(model_MCMC_WI_Sol)[grepl('pop_label', colnames(model_MCMC_WI_Sol))]
  breeding_values <- data.frame(model_MCMC_WI_Sol[, col_subset])
  new_colnames <- gsub(".*pop_label\\.", "", colnames(breeding_values))
  colnames(breeding_values) <- new_colnames
  
  # Filter genotype matrix M to match the breeding values' column names
  M <- M[rownames(M) %in% new_colnames, ]
  
  ## Backsolving process ##
  
  # Calculate the transpose of M
  M_transpose <- t(M)
  
  # Calculate the product MM'
  MM_prime <- M %*% M_transpose
  
  # Calculate the inverse of MM'
  MM_prime_inverse <- solve(MM_prime)
  
  # Convert breeding values to matrix
  breeding_values_matrice <- as.matrix(breeding_values)
  
  # Calculate SNP effects
  SNPs_effects <- M_transpose %*% MM_prime_inverse %*% t(breeding_values_matrice)
  
  # Write the SNP effects to a CSV file
  write.csv(data.frame(cbind(snp_positions, SNPs_effects)),
            file = gsub('LLL', output, gsub('XXX', p, 'SNPs_effects_XXX_LLL.csv')), quote = FALSE)
  
  # Write the breeding values to a CSV file
  write.csv(breeding_values,
            file = gsub('LLL', output, gsub('XXX', p, 'BreedingValues_XXX_LLL.csv')), quote = FALSE)
}
