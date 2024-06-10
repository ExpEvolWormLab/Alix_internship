### Isotypes ###
# Script to create isotype group (Strain pairs with concordance > 0.99985)
# First and Second part should be split from others, and just produce the matrix in a .tsv file

# To change
# working directory
setwd("~/Documents/Worms/VariantCalling")
# hard filtering file after imputation
file <- "AllLines.hard_filter.vcf.gz"
# Cutoff
T <- 0.99985

library(vcfR)
library(ggplot2)

##### FIRST PART #####
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

nbr_na <- data.frame(NA_values=colSums(is.na(convert_genotype))/nrow(convert_genotype))

write.csv(nbr_na,'Na_distribution.tsv',sep = '\t',row.names = FALSE, quote = FALSE)