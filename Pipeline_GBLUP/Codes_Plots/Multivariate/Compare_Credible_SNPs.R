### MULTI-TRAIT
## Output : 
#    Shared_credible_SNPs*tsv : file with two columns Comparison Nbr_SNPs, compare each traits and store the number of credible SNPs they share 
#    Shared_credible_SNPs*pdf : plot the previous files
#    SNP_Count : file which store for SNPs the number of group it's associated to 
#    Repartition_SNPs_along_K : for each SNPs plot it along chromosome in function of the number of groups they are associated with

### Script which compare the credible SNP for each trait
output <- 'VanRaden_A6_NGM_0.99'
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Pruned")

# Load necessary libraries
library(tidyr)
library(ggplot2)

# Load the CSV files matching the pattern
motif <- paste0("Summary_.*_", output, ".csv")
motif1 <- paste0("Summary_(.*)_", output, ".csv")
files <- list.files(pattern = motif)
file_names <- sub(motif1, "\\1", files)
dfs <- lapply(files, read.csv)
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results")

# Filter for rows marked as "Credible"
credibles <- lapply(dfs, function(df) df[df$credible == "Credible",]$X)

# Function to compare the names of credible rows between groups of data frames
compare_lists <- function(lists) {
  if (length(lists) == 1) {
    return(lists[[1]])
  } else {
    result <- lists[[1]]
    for (i in 2:length(lists)) {
      result <- intersect(result, lists[[i]])
    }
    return(result)
  }
}

# Generate all possible combinations
library(combinat)

# Store results for individual comparisons
results_solo <- list()
results_solo[[paste0("combinations_", 1)]] <- lapply(seq_along(dfs), function(idx) {
  list(
    combination = file_names[idx],
    intersection = length(credibles[[idx]])
  )
})

# Store results for combined comparisons
results <- list()
for (i in 2:length(credibles)) {
  comb <- combn(seq_along(credibles), i, simplify = FALSE)
  comb_names <- combn(file_names, i, simplify = FALSE)
  results[[paste0("combinations_", i)]] <- lapply(seq_along(comb), function(idx) {
    list(
      combination = comb_names[[idx]],
      intersection = compare_lists(credibles[comb[[idx]]])
    )
  })
}

# Write results to a file
output_file <- file(gsub('XXX', output, 'Shared_credible_SNPs_XXX.tsv'), open = "wt")

# Write header to the file
cat("Comparison\tNbr_SNPs\n", file = output_file)

# Write solo results to the file
for (i in seq_along(results_solo)) {
  for (k in seq_along(results_solo[[i]])) {
    cat(paste(results_solo[[i]][[k]]$combination, collapse = ","), "\t", file = output_file)
    cat(results_solo[[i]][[k]]$intersection, "\n", file = output_file)
  }
  cat("\n", file = output_file)
}

# Write combined results to the file
for (i in seq_along(results)) {
  print(i)
  if (i == 6) {
    cat(paste(results[[i]]$combination, collapse = ","), "\t", file = output_file)
    cat(length(results[[i]]$intersection), "\n", file = output_file)
    break
  }
  for (k in seq_along(results[[i]])) {
    cat(paste(results[[i]][[k]]$combination, collapse = ","), "\t", file = output_file)
    cat(length(results[[i]][[k]]$intersection), "\n", file = output_file)
  }
  if (i < length(results)) {
    cat("\n", file = output_file)
  }
}

# Close the output file
close(output_file)

# Read the output file into a data frame
data <- read.delim(gsub('XXX', output, 'Shared_credible_SNPs_XXX.tsv'))

# Add a column indicating the number of elements compared in each combination
data$Num_Comparisons <- sapply(strsplit(as.character(data$Comparison), ","), length)

# Plot the number of credible SNPs shared between different traits
ggplot(data, aes(x = Comparison, y = Nbr_SNPs, fill = Comparison)) +
  geom_bar(stat = "identity") +
  labs(x = "Traits", y = "Number of SNPs", title = "Number of credible SNPs Shared between Different Trait") +
  theme_minimal() +
  scale_x_discrete(labels = function(x) gsub(",", "\n", x)) +
  facet_wrap(~ Num_Comparisons, ncol = 2, scales = "free") +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))

# Save the plot to a PDF file
ggsave(gsub('XXX', output, 'Shared_credible_SNPs_XXX_plot.pdf'))

# Write a file containing the SNP ID and the number of traits each credible SNP is associated with
# Create a list to store the count of each SNP
snp_count <- list()

# Loop through each data frame to count credible SNPs
for (i in 1:length(dfs)) {
  snps <- dfs[[i]][dfs[[i]]$credible == "Credible", "X"]
  for (snp in snps) {
    if (is.null(snp_count[[snp]])) {
      snp_count[[snp]] <- 1
    } else {
      snp_count[[snp]] <- snp_count[[snp]] + 1
    }
  }
}

# Create a data frame with SNP and their counts
snp_df <- data.frame(
  SNP = names(snp_count),
  Count = unlist(snp_count)
)

# Split the SNP column into separate CHROM and POS columns
snp_df <- separate(snp_df, SNP, into = c("CHROM", "POS"), sep = "_", remove = FALSE)

# Write the data frame to a CSV file
write.csv(snp_df, gsub('XXX', output, "SNP_Count_XXX.csv"), quote = FALSE, row.names = FALSE)

# Find SNPs that are not credible and set their count to 0
all_snps <- dfs[[1]][,'X']
not_credible <- data.frame(SNP = all_snps[!(all_snps %in% snp_df$SNP)], Count = 0)
not_credible <- separate(not_credible, SNP, into = c("CHROM", "POS"), sep = "_", remove = FALSE)

# Combine credible and not credible SNPs into one data frame
SNP_df <- data.frame(rbind(not_credible, snp_df))

# Manhattan plot for the distribution of SNPs along chromosomes
# Convert the POS variable to numeric
SNP_df$POS <- as.numeric(SNP_df$POS)

# Create the plot using ggplot2
repartition <- ggplot(SNP_df, aes(x = POS, y = Count, color = factor(Count))) +
  geom_point(size = 0.5) +
  theme_minimal() +
  labs(title = "Repartition SNPs along Chromosomes",
       x = "SNPs",
       y = "Number of trait associated") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        axis.ticks.x = element_blank()) +
  facet_wrap(~ CHROM, ncol = 2, scales = "free_x") +
  scale_x_continuous(breaks = seq(0, max(SNP_df$POS), by = 1000000)) +
  theme(legend.position = "none")

# Print the plot
print(repartition)

# Save the plot to a PDF file
ggsave(gsub('XXX', output, 'Repartition_SNPs_along_K_XXX.pdf'))
