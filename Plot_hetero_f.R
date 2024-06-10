### Script which compute and plot heterozygosity along chromosome for each population

### User's variables
# Path where data are
#setwd("~/Documents/Worms/VariantCalling/Hetero")
# File to analyze (multi-sampling vcf)
#file <- "II.final.vcf.gz"
file <- commandArgs(trailingOnly = TRUE)
if (length(file) == 0) {
  stop("No VCF file provided. Please specify the VCF file to analyze.")
}
# Chromosomes to plot
chromosomes <- c('I', 'II', 'III', 'X')
#chromosomes <- c('I', 'II', 'III', 'IV', 'V', 'X')
# Population to plot
populations <- c("A6", "CA1", "CA2", "CA3", "EEV", "GA1", "GA2", "GA4", "GM1", "GM3", "GT1", "GT2", "LR1", "LR3", "SMR2", "SMR4")
populations <- c("A6", "CA", "EEV", "GA",  "GM", "GT","LR",  "SMR")
#populations <- c(c('EEV1401','SMR24L102'))
populations <- c('*')
# Groups to make
group1 <- c("A6", "EEV")
group2 <- c("CA1", "CA2", "CA3")
group3 <- c("GA1", "GA2", "GA4")
group4 <- c("GM1", "GM3", "GT1", "GT2")
group5 <- c("LR1", "LR3", "SMR2", "SMR4")
groups <- list(group1, group2, group3, group4, group5)
list_group1 <- list()
list_group2 <- list()
list_group3 <- list()
list_group4 <- list()
list_group5 <- list()
list_groups <- list(list_group1, list_group2, list_group3, list_group4, list_group5)
#name_to_plot <- c('EEV1401','SMR24L102')
# Define number of bins per chromosome
size_bin <- 500

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

start_time <- Sys.time()
last_bin <- 0
Results <- data.frame(Chromosome = character(0), Bin_Number = numeric(0), Start = character(0), End = character(0), Number_of_heterozygotes = numeric(0))

# Loop through each chromosome
for (pop in populations) {
  for (i in chromosomes) {
    print(paste0("Start of chromosome ", i))
    
    # Extract df of this chromosome
    df <- snp_genotypes[snp_genotypes$CHROM == i, ]
    
    # Remove CHROM column
    df$CHROM <- NULL
    
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
  g <- ggplot() +
    geom_point(data = Results, aes(x = Bin_Number, y = Number_of_heterozygotes, color = Chromosome)) +
    ggtitle(paste0("Heterozygosity along chromosome for ", pop)) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  if (pop %in% group1) {list_group1[[length(list_group1) + 1]] <- g}
  if (pop %in% group2) {list_group2[[length(list_group2) + 1]] <- g}
  if (pop %in% group3) {list_group3[[length(list_group3) + 1]] <- g}
  if (pop %in% group4) {list_group4[[length(list_group4) + 1]] <- g}
  if (pop %in% group5) {list_group5[[length(list_group5) + 1]] <- g}
  
  # Save plot as PDF
  pdf_filename <- paste0("Heterozygosity_", gsub(" ", "_", pop), ".pdf")
  ggsave(pdf_filename, g)
  
  # Display plot
  print(g)
  dev.off()
  
  # Clear memory
  rm(g)
  rm(results)
  write.csv(Results, "Heterozygosity.csv", row.names = FALSE, quote = FALSE)
  # Clear Results for the next population
  Results <- data.frame(Chromosome = character(0), Bin_Number = numeric(0), Start = character(0), End = character(0), Number_of_heterozygotes = numeric(0))
  
  # Break if you only want to generate plots for one population
  #break
}



## Plot the plot per groups
for (i in seq_along(groups)) {
  pdf(file = paste0("Plots_for_", paste(groups[[i]], collapse = "_"), ".pdf"))
  grid.arrange(grobs = list_groups[[i]], ncol = 2)
  dev.off()
}

# Record the end time
end_time <- Sys.time()

# Calculate the time taken
time_taken <- end_time - start_time

# Print the time taken
print(time_taken)


Results[is.na(Results)] <- 0

df <- read.delim("20220216_c_elegans_divergent_regions_strain.bed", header = FALSE)
df$V4 <- NULL
df <- unique(df)

# Function to merge overlapping intervals
merge_intervals <- function(df) {
  result <- list()
  current_interval <- c(df$V2[1], df$V3[1])
  
  for (i in 2:nrow(df)) {
    if (df$V2[i] <= current_interval[2]) {
      current_interval[2] <- max(current_interval[2], df$V3[i])
    } else {
      result <- rbind(result, current_interval)
      current_interval <- c(df$V2[i], df$V3[i])
    }
  }
  
  result <- rbind(result, current_interval)
  return(as.data.frame(result))
}

# Apply function to merge intervals
chrom <- unique(df$V1)

HD_df <- data.frame(Chromosome = character(0), Start = integer(0), End = integer(0))
for (i in chrom) {
  df_chrom <- df[df$V1 == i, ]
  res <- merge_intervals(df_chrom)
  colnames(res) <- c("Start", "End")
  res$Chromosome <- rep(i, nrow(res))
  HD_df <- rbind(HD_df, res)
}

row.names(HD_df) <- NULL

# Initialize HD_r vector with the length of Results$Chromosome
HD_r <- integer(length(Results$Chromosome))
# Initialize nbr variable to keep track of the number of rows processed
nbr <- 0

# Loop over each chromosome in the chrom vector
for (c in chrom) {
  # Subset Results and HD_df data frames based on the current chromosome
  Results_chrom <- Results[Results$Chromosome == c, ]
  HD_df_chrom <- HD_df[HD_df$Chromosome == c, ]
  
  # Loop over each row of HD_df_chrom to iterate through the regions
  for (j in 1:nrow(HD_df_chrom)) {
    # Extract start and end coordinates of the current region
    start <- HD_df_chrom$Start[j]
    end <- HD_df_chrom$End[j]
    
    # Loop over each row of Results_chrom to check for overlap with the current region
    for (i in 1:nrow(Results_chrom)) {
      # Check if there is overlap between the current region and the current result
      if (start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
          start <= Results_chrom$Start[i] && end >= Results_chrom$End[i] || 
          start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
          start >= Results_chrom$Start[i] && end <= Results_chrom$End[i]) {
        # If there is overlap, mark the corresponding entry in HD_r as -1
        HD_r[i + nbr] <- -0.005
      } else {
        # If there is no overlap, mark the corresponding entry in HD_r as NA
        HD_r[i + nbr] <- NA
      }
    }
  }
  
  # Update nbr to keep track of the total number of rows processed
  nbr <- nbr + nrow(Results_chrom)
}
Results$HD <- HD_r

g <- ggplot() +
  geom_point(data = Results, aes(x = Bin_Number, y = Number_of_heterozygotes, color = Chromosome)) +
  geom_point(data = Results, aes(x = Bin_Number, y = HD)) +
  ggtitle(paste0("Heterozygosity along chromosome for ", pop)) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#print(g)

Results[is.na(Results)] <- 0

# Sort dataframe
Results_sorted <- Results[order(-Results$Number_of_heterozygotes), ]
rownames(Results_sorted) <- c(1:length(Results$Chromosome))
# Get the 5% most important values
nbr <- 5 * length(Results$Chromosome) / 100
Results_sorted <- Results_sorted[1:nbr, ]
#Results_sorted[is.na(Results_sorted)] <- 0
prop_hetero <- sum(-Results_sorted$HD) / nbr

### RANDOMISATION IN ORDER TO SEE IF SIGNIFICANT OR NOT
randomisation <- function(df, nbr) {
  df_random <- as.data.frame(df$HD)
  colnames(df_random) <- "HD"
  df_random$Number_of_heterozygotes <- sample(df$Number_of_heterozygotes)
  df_sorted <- df_random[order(-df_random$Number_of_heterozygotes), ]
  rownames(df_sorted) <- c(1:length(df_sorted$Number_of_heterozygotes))
  df_sorted <- df_sorted[1:nbr, ]
  #df_sorted[is.na(df_sorted)] <- 0
  prop_hetero <- sum(-df_sorted$HD) / nbr
  return(prop_hetero)
}

# Set the number of repetitions
num_repetitions <- 50000

# Apply the randomisation function 1000 times and store the results in a vector
distribution <- replicate(num_repetitions, randomisation(Results, nbr))

hist <- ggplot(NULL, aes(x = distribution)) +
  geom_histogram(color = "black", fill = "white") +
  geom_vline(xintercept = prop_hetero, color = "red") +
  theme_bw()

pdf_filename <- paste0("Randomisation_", gsub(" ", "_", pop), ".pdf")
ggsave(pdf_filename, hist)

print(hist)
