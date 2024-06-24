### Script to create isotype group (Strain pairs with concordance > 0.99985)
### First and Second part should be split from others, and just produce the matrix in a .tsv file

# Change the working directory
setwd("~/Documents/Worms/VariantCalling/Isotypes")

# Hard filtering file after imputation
file <- "I.final.hard_filter.vcf.gz"

# Proportion of NaN cutoff
T <- 0.982

# Load necessary libraries
library(vcfR)
library(ggplot2)

### Read vcf file
### Convert genotype in number (-1 for Na, 0 or 2)

start_time <- Sys.time()

# Read vcf file
vcf <- read.vcfR(file)

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]

# Extract genotype information
genotype_info <- vcf@gt

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

##### Filter lines with too many NA #####

# Read the distribution of NA values
nbr_na <- read.csv('Na_distribution.tsv')
rownames(nbr_na) <- nbr_na$X
nbr_na$X <- NULL

# Plot the histogram of missing values distribution for all strains
g <- ggplot(nbr_na, aes(x = NA_values)) + 
  geom_histogram(fill = 'cyan', color = 'black', binwidth = 0.007) + # each bin should have a width of 0.007 unit
  geom_vline(xintercept = T, color = 'red') +
  labs(title = "Missing values distribution", x = "Values", y = "Number") +
  theme_bw()

# Save the plot
ggsave('histo_NA.pdf', g)

# Filter strains based on the NA values threshold
to_keep <- rownames(nbr_na[nbr_na$NA_values < T,, drop = FALSE]) # drop prevents removal of rownames

# Reduce the genotype matrix to only include the filtered strains
reduced_convert_genotype <- convert_genotype[, to_keep]

### Compute matrix which stores how many SNPs are similar between two lines

# Number of SNPs
N_snp <- nrow(reduced_convert_genotype)

# Number of individuals (strains)
N_ind <- ncol(reduced_convert_genotype)

# Initialize the similarity matrix with zeros
matrix <- matrix(0, nrow = ncol(reduced_convert_genotype), ncol = ncol(reduced_convert_genotype))

# Compare the genotypes pairwise
for (i in 1:N_ind) {
  progress_percentage <- (i / N_ind) * 100
  cat(sprintf("\rProgress: %.2f%%", progress_percentage))
  for (j in (i+1):(N_ind-1)) {
    subset <- na.omit(reduced_convert_genotype[, c(i, j)])
    matrix[i, j] <- sum(subset[, 1] == subset[, 2]) / nrow(subset)
    matrix[j, i] <- matrix[i, j] 
  }
}

# Set row and column names of the matrix
colnames(matrix) <- colnames(reduced_convert_genotype)
rownames(matrix) <- colnames(reduced_convert_genotype)

# Record the end time for performance measurement
end_time <- Sys.time()

# Calculate and print the time taken for the computation
time_taken <- end_time - start_time
print(time_taken)

# Write the similarity matrix to a .tsv file
write.table(matrix, paste0("matrix_concordance.by_hand.",T,".tsv"), quote = FALSE, sep = '\t')
