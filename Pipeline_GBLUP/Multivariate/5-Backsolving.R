## Script to backsolve multivariate analysis ##

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]  # Path to the Sol file
file1 <- args[2]  # Path to the M file
file2 <- args[3]  # Path to the SNP positions file

# Name of the output directory
output <- args[4]

# Load necessary library
library(data.table)

# Load model solutions and M matrix
model_MCMC_WI_Sol <- read.csv(file)
M <- as.data.frame(fread(file1))
rownames(M) <- M[, 1]  # Set row names
M[, 1] <- NULL  # Remove the first column
M <- as.matrix(M)  # Convert to matrix

# Load SNP positions
snp_positions <- read.csv(file2)

# Define the vector of traits
vect_P_traits <- c("SF",
                   "SB",
                   "FS",
                   "FB",
                   "BS",
                   "BF")

# Loop through each trait to perform backsolving
for(p in vect_P_traits) {
  # Get the vector of breeding values for each trait
  col_subset <- colnames(model_MCMC_WI_Sol)[grepl(paste0(p, '.pop_label'), colnames(model_MCMC_WI_Sol))]
  breeding_values <- data.frame(model_MCMC_WI_Sol[, col_subset])
  
  # Update column names to remove prefix
  new_colnames <- gsub(".*pop_label\\.", "", colnames(breeding_values))
  colnames(breeding_values) <- new_colnames
  
  # Filter M matrix to include only relevant columns
  M <- M[rownames(M) %in% new_colnames, ]
  
  # Backsolving
  # Compute the transpose of M
  M_transpose <- t(M)
  
  # Compute MM'
  MM_prime <- M %*% M_transpose
  
  # Compute the inverse of MM'
  MM_prime_inverse <- solve(MM_prime)
  
  # Convert breeding values to matrix
  breeding_values_matrice <- as.matrix(breeding_values)
  
  # Compute SNP effects
  SNPs_effects <- M_transpose %*% MM_prime_inverse %*% t(breeding_values_matrice)
  
  # Write SNP effects to CSV
  write.csv(data.frame(cbind(snp_positions, SNPs_effects)), 
            file = gsub('LLL', output, gsub('XXX', p, 'SNPs_effects_XXX_LLL.csv')), 
            quote = FALSE)
  
  # Write breeding values to CSV
  write.csv(breeding_values, 
            file = gsub('LLL', output, gsub('XXX', p, 'BreedingValues_XXX_LLL.csv')), 
            quote = FALSE)
}
