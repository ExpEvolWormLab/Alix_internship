### Script to summarize SNPs effects
## Input : 1. SNPs_effects*csv - output of Pipeline_GBLUP
##         2. Name of output
## Output : Summary*csv - Name of SNP (CHROM_POS), CHROM, Median (of SNPs effect a posteriori distribution), lower (of credibility interval), upper, Credibility

# Load necessary libraries
library(HDInterval)  # For computing highest density intervals
library(data.table)  # For fast data reading and manipulation
library(ggplot2)  # For plotting
library(matrixStats)  # For efficient matrix operations

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]  # The input file containing SNP effects
output <- args[2]  # The output file name

# Read the input file into a data frame
Df <- as.data.frame(fread(file))
print('Data frame loaded')

# Set row names of the data frame to the first column and remove the first column
row.names(Df) <- Df$V1
Df$V1 <- NULL

# Extract chromosome information and remove the CHROM and POS columns from the data frame
chrom <- Df$CHROM
Df$CHROM <- NULL
Df$POS <- NULL

# Compute the 98% highest density intervals (HDI) for the SNP effects
HDI <- hdi(t(Df), 0.98)
print('Credibility interval computed')

# Function to check if zero is within the interval
is_zero_in_interval <- function(lower, upper) {
  return(lower <= 0 & upper >= 0)
}

# Identify credible SNPs where zero is not in the interval
sig <- colnames(HDI)[!apply(HDI, 2, function(x) is_zero_in_interval(x["lower"], x["upper"]))]
print('Credible SNPs found')

# Add columns for the median and the credibility intervals
Df$median <- rowMedians(as.matrix(Df))
Df$lower <- HDI['lower',]
Df$upper <- HDI['upper',]

# Add a column to indicate credible SNPs
Df$credible <- ifelse(rownames(Df) %in% sig, "Credible", "Not Credible")

# Add back the chromosome information
Df$CHROM <- as.factor(chrom)

# Write the summary to a CSV file
write.csv(Df[, c('CHROM', 'median', 'lower', 'upper', 'credible')], gsub('XXX', output, 'Summary_XXX.csv'))
