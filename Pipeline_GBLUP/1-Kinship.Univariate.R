## Script to compute the kinship matrix in three different ways (Hoffman, VanRaden, Noia)

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]  # VCF file
output <- args[2]  # Output file base name

# Load necessary packages
library(vcfR)
library(ggplot2)
library(rutilstimflutre)

##### First step: kinship matrix ######

# Read VCF file
vcf <- read.vcfR(file)

# Extract genotype information
genotype_info <- vcf@gt

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]

# Remove FORMAT column if it exists
if ("FORMAT" %in% colnames(genotype_info)) {
  genotype_info <- subset(genotype_info, select = -c(FORMAT))
}

# Function to convert genotypes to numeric values
get_genotype <- function(genotype) {
  ifelse(substr(genotype, 1, 3) %in% c("0/0", "0|0"), 0,
         ifelse(substr(genotype, 1, 3) %in% c("1/1", "1|1"), 2,
                ifelse(substr(genotype, 1, 3) %in% c("./.", ".|."), NaN, NaN)))
}

# Apply the conversion function to genotype data
convert_genotype <- apply(genotype_info, 2, get_genotype)

# Count NA values per SNP
nbr_na <- data.frame(NA_values = rowSums(is.na(convert_genotype)))

# Plot the distribution of NA values per SNP
g <- ggplot(nbr_na, aes(x = NA_values)) +
  geom_histogram() +
  theme_bw() +
  xlim(c(1, 20)) +
  ggtitle('NA distribution per SNPs')

# Save the NA distribution plot
ggsave(filename = gsub('XXX', output, 'XXX_NA_distribution.pdf'), plot = g)

# Remove rows with NA values from SNP positions and genotypes
snp_positions <- as.data.frame(snp_positions[complete.cases(convert_genotype), ])
convert_genotype <- convert_genotype[complete.cases(convert_genotype), ]

# Transpose the genotype matrix
M <- t(convert_genotype)

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

# Clean row and column names of the genotype matrix
rownames(M) <- clean_colnames(rownames(M))
rownames(M) <- gsub('CeMee', '', rownames(M))
rownames(M) <- gsub('_sorted', '', rownames(M))
colnames(M) <- paste0(snp_positions$CHROM, '_', snp_positions$POS)

## VanRaden 2008 kinship matrix
VR <- estimGenRel(
  M,
  relationships = "additive",
  method = "vanraden1"
)

## Noia (Vitezica et al, 2017) kinship matrix
noia <- estimGenRel(
  M,
  relationships = "additive",
  method = "noia"
)

## Hoffman method for calculating genetic similarity
calcular_similarite_genetique <- function(gmat) {
  m <- nrow(gmat)  # Number of markers
  gmat <- scale(t(gmat))  # Scale and transpose the genotype matrix
  gmat[is.na(gmat)] <- 0  # Replace NA values with 0
  X <- tcrossprod(gmat) / m  # Calculate the cross-product matrix
  svdx <- svd(X)  # Perform singular value decomposition
  D <- diag(sqrt(svdx$d))  # Calculate diagonal matrix of singular values
  US <- svdx$u %*% D  # Calculate US matrix
  K <- tcrossprod(US)  # Calculate the kinship matrix
  K <- K / mean(diag(K))  # Normalize the kinship matrix
  rownames(K) <- colnames(K) <- rownames(gmat)  # Set row and column names
  return(K)
}

# Calculate kinship matrix using the Hoffman method
A <- calcular_similarite_genetique(convert_genotype)

# Save the results to CSV files
write.csv(snp_positions, gsub('XXX', output, 'XXX_snp_positions.csv'), row.names = FALSE, quote = FALSE)  # SNP positions
write.csv(M, gsub('XXX', output, 'XXX_t_convert_genotype.csv'), quote = FALSE)  # Transposed genotype matrix
write.csv(VR, gsub('XXX', output, 'Kinship_matrix_VanRaden_XXX.csv'), row.names = FALSE, quote = FALSE)  # VanRaden kinship matrix
write.csv(A, gsub('XXX', output, 'Kinship_matrix_Luke_XXX.csv'), row.names = FALSE, quote = FALSE)  # Hoffman kinship matrix
write.csv(noia, gsub('XXX', output, 'Kinship_matrix_noia_XXX.csv'), row.names = FALSE, quote = FALSE)  # Noia kinship matrix
