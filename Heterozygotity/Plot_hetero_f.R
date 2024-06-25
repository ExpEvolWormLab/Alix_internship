### Script which compute heterozygotity along chromosome and store it in a .tsv file

# Chromosomes to plot
chromosomes <- c('I', 'II', 'III', 'IV', 'V', 'X')
# Population to plot
populations <- c("A6", "CA", "EEV", "GA",  "GM", "GT","LR",  "SMR")
# Define number of SNPs per bin
size_bin <- 500

# Get the VCF file from command line arguments
file <- commandArgs(trailingOnly = TRUE)
if (length(file) == 0) {
  stop("No VCF file provided. Please specify the VCF file to analyze.")
}

# Load necessary libraries
library(vcfR)
library(ggplot2)
library(gridExtra)

# Read VCF file
vcf <- read.vcfR(file)

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]

# Extract genotype information
genotype_info <- vcf@gt

# Define a function to determine homozygosity/heterozygosity
get_genotype <- function(genotype) {
  ifelse(grepl("is_het", genotype), 1, 0)
}

# Determine genotype for each sample
sample_genotypes <- apply(genotype_info, 2, get_genotype)

# Remove FORMAT column if it exists
if ("FORMAT" %in% colnames(sample_genotypes)) {
  sample_genotypes <- subset(sample_genotypes, select = -c(FORMAT))
}

# Combine SNP positions and genotype information
snp_genotypes <- as.data.frame(cbind(snp_positions, sample_genotypes))
snp_genotypes$POS <- as.integer(snp_genotypes$POS)

# Function to compute heterozygosity in bins for a population
compute_hist_pop <- function(df, list_bin, iter, pop, last_bin, chrom) {
  # Initialize vectors to store results
  length_het <- numeric(iter - 1)
  list_bin_start <- numeric(iter - 1)
  list_bin_end <- numeric(iter - 1)
  list_nbr_bin <- numeric(iter - 1)
  list_chrom <- rep(chrom, iter - 1)  # Repeat chromosome name for each bin
  
  # Loop over bins
  for (i in seq_len(iter - 1)) {
    progress_percentage <- (i / iter) * 100
    cat(sprintf("\rProgress: %.2f%%", progress_percentage))
    
    # Save bin number and chromosome for this iteration
    list_nbr_bin[i] <- i + last_bin
    
    # Select subset of dataframe within current bin interval
    df1 <- df[(list_bin[i] * size_bin):(list_bin[i + 1] * size_bin), ]
    # Remove POS column and store start and end
    list_bin_start[i] <- df1[1, ]$POS
    list_bin_end[i] <- df1[nrow(df1), ]$POS
    df1 <- df1[, -which(names(df1) == "POS")]
    # Compute mean heterozygosity for each individual in the bin
    heterozygous_counts <- rowMeans(df1 == 1)
    # Mean heterozygosity across individuals in the bin
    length_het[i] <- mean(heterozygous_counts)
  }
  
  # Combine results into a data frame
  results <- data.frame(Chromosome = list_chrom, Bin_Number = list_nbr_bin, Start = list_bin_start, End = list_bin_end, Number_of_heterozygotes = length_het)
  return(results)
}

# Start time for performance measurement
start_time <- Sys.time()
last_bin <- 0
Results <- data.frame(Chromosome = character(0), Bin_Number = numeric(0), Start = character(0), End = character(0), Number_of_heterozygotes = numeric(0))

# Loop through each population and chromosome
for (pop in populations) {
  for (i in chromosomes) {
    print(paste0("Start of chromosome ", i))
    
    # Extract dataframe of this chromosome
    df <- snp_genotypes[snp_genotypes$CHROM == i, ]
    
    # Remove CHROM column
    df$CHROM <- NULL
    
    # Define bins and iterations
    bin <- size_bin
    nbr_snp <- nrow(df)
    list_bins <- seq(0, nbr_snp / bin)
    iter <- length(list_bins)
    
    # Select columns we want to plot
    columns <- colnames(df)[grepl(pop, colnames(df))]
    df1 <- df[, c("POS", columns)]
    results <- compute_hist_pop(df1, list_bins, iter, pop, last_bin, i)
    last_bin <- as.numeric(results[nrow(results), "Bin_Number"])
    Results <- rbind(Results, results)
    print(paste0("End of chromosome ", i))
  }

  # Clear memory
  rm(g)
  rm(results)
  write.csv(Results, gsub('XXX',pop,"Heterozygosity_XXX.csv"), row.names = FALSE, quote = FALSE)
  # Clear Results for the next population
  Results <- data.frame(Chromosome = character(0), Bin_Number = numeric(0), Start = character(0), End = character(0), Number_of_heterozygotes = numeric(0))
  
  # Break if you only want to generate plots for one population
  #break
}

# Record the end time
end_time <- Sys.time()

# Calculate the time taken
time_taken <- end_time - start_time

# Print the time taken
print(time_taken)
