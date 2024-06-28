## Script to backsolve with pruned matrix for UNIVARIATE ##
# Define the output, vcf file, and condition
output <- 'Uni_VanRaden_A6_NaCl_0.99'
output_dir <- "~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned"
file1 <- 'pruned.0.99.vcf.gz'
condition <- 'NaCl'

# Set the working directory to where the pruned VCF file is located
setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")

# Define the traits
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")

# Load necessary libraries
library(data.table)
library(vcfR)

#

# Read the VCF file
vcf <- read.vcfR(file1)

# Extract genotype information
genotype_info <- vcf@gt

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]

# Remove FORMAT column if it exists
if ("FORMAT" %in% colnames(genotype_info)) {
  genotype_info <- subset(genotype_info, select = -c(FORMAT))
}

# Function to convert genotypes to numerical values
get_genotype <- function(genotype) {
  ifelse(substr(genotype, 1, 3) %in% c("0/0", "0|0"), 0,
         ifelse(substr(genotype, 1, 3) %in% c("1/1", "1|1"), 2,
                ifelse(substr(genotype, 1, 3) %in% c("./.", ".|."), NaN, NaN)))
}

# Apply the genotype conversion function to the genotype information
convert_genotype <- apply(genotype_info, 2, get_genotype)

# Filter out rows with missing genotype information
snp_positions <- as.data.frame(snp_positions[complete.cases(convert_genotype), ])
convert_genotype <- convert_genotype[complete.cases(convert_genotype), ]

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

# Clean the column and row names of the genotype matrix
colnames(convert_genotype) <- clean_colnames(colnames(convert_genotype))
M <- as.matrix(t(convert_genotype))
rownames(M) <- gsub('CeMee', '', rownames(M))
rownames(M) <- gsub('_sorted', '', rownames(M))

# Loop through each trait and perform the backsolving
for (p in vect_P_traits) {
  # Set the working directory to where the model solutions are stored
  setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
  
  # Define the file name for the model solution based on the trait and condition
  file <- gsub('condition', condition, gsub('TTT', p, 'VanRaden_A6_condition_TTT_MCMCmodel_Sol.csv'))
  
  # Read the model solution file
  model_MCMC_WI_Sol <- read.csv(file)

  # Extract the breeding values for each trait
  col_subset <- colnames(model_MCMC_WI_Sol)[grepl('pop_label', colnames(model_MCMC_WI_Sol))]
  breeding_values <- data.frame(model_MCMC_WI_Sol[, col_subset])
  new_colnames <- gsub(".*pop_label\\.", "", colnames(breeding_values))
  colnames(breeding_values) <- new_colnames
  M <- M[rownames(M) %in% new_colnames, ]
  
  # Backsolving process
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
  
  # Set the working directory to where the results should be saved
  setwd(output_dir)
  
  # Save the SNP effects and breeding values to CSV files
  write.csv(data.frame(cbind(snp_positions, SNPs_effects)), file = gsub('LLL', output, gsub('XXX', p, 'SNPs_effects_XXX_LLL.csv')), quote = FALSE)
  write.csv(breeding_values, file = gsub('LLL', output, gsub('XXX', p, 'BreedingValues_XXX_LLL.csv')), quote = FALSE)
}
