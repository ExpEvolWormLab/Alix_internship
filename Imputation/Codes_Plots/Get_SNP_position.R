## Script to get the position of SNPs present in a set of VCF files with similar names

# Set working directory
setwd("~/Documents/Worms/Plot_GATK/Pop")

# Define the pattern to match VCF files
motif <- '(.*)_final.vcf.gz'

# Load necessary library
library(vcfR)

# List all files in the directory matching the pattern
files <- list.files(pattern = motif)

# Read each VCF file into a list of vcfR objects
dfs <- lapply(files, read.vcfR)

# Extract the base names of the files (without the suffix)
files_name <- sub(motif, "\\1", files)

# Loop through each VCF object and extract SNP positions
for(i in seq_along(dfs)){
  vcf <- dfs[[i]]  # Get the i-th vcfR object
  
  # Extract SNP positions (CHROM and POS columns from the fixed section of the VCF)
  snp_positions <- vcf@fix[, c("CHROM", "POS")]
  
  # Write the SNP positions to a CSV file
  write.csv(snp_positions, paste0(files_name[i], '_snp_positions.csv'), row.names = FALSE, quote = FALSE)
}
