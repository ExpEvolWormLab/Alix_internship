### Script to plot phenotype and genotype data as in Stephens and al. (2013)
## Need to specify a specific position in genome

# Set the working directory
setwd("~/Documents/Worms/GBLUP")

# Define file paths and output names
matrix <- 'Inverted_kinship_matrix_VanRaden_A6.csv'
output <- 'A6_NaCl_NGM'
file <- 'pruned.0.99.vcf.gz'
chrom <- 'II'
pos <- '3373318'


# Load necessary libraries
library(vcfR)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Read the VCF file
vcf <- read.vcfR(file)

# Extract genotype information
genotype_info <- vcf@gt

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]

# Remove FORMAT column if it exists
if ("FORMAT" %in% colnames(genotype_info)) {
  genotype_info <- subset(genotype_info, select = -c(FORMAT))
}

# Function to convert genotypes
get_genotype <- function(genotype) {
  ifelse(substr(genotype, 1, 3) %in% c("0/0", "0|0"), 0,
         ifelse(substr(genotype, 1, 3) %in% c("1/1", "1|1"), 2,
                ifelse(substr(genotype, 1, 3) %in% c("./.", ".|."), NaN, NaN)))
}

# Function to clean column names
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

# Convert genotype information
convert_genotype <- apply(genotype_info, 2, get_genotype)
# Clean the column and row names
colnames(convert_genotype) <- clean_colnames(colnames(convert_genotype))
colnames(convert_genotype) <- gsub('CeMee', '', colnames(convert_genotype))
colnames(convert_genotype) <- gsub('_sorted', '', colnames(convert_genotype))
convert_genotype <- data.frame(cbind(snp_positions, convert_genotype))

# Select a specific SNP position
subset_genotype <- convert_genotype[convert_genotype$CHROM == chrom & convert_genotype$POS == pos,]
subset_genotype$CHROM <- NULL
subset_genotype$POS <- NULL
subset_genotype <- as.data.frame(t(subset_genotype))
subset_genotype$pop_label <- rownames(subset_genotype)
colnames(subset_genotype) <- c('Genotype', 'pop_label')

# Load and process kinship matrix
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP")
Ginv <- as.matrix(read.csv(matrix))
# Corrected names
colnames(Ginv) <- gsub('CeMee', '', colnames(Ginv))
colnames(Ginv) <- gsub('_sorted', '', colnames(Ginv))
rownames(Ginv) <- colnames(Ginv)

# Load phenotype data
pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
pheno <- pheno[, c('pop_label', 'temperature', 'rel_humidity', "session_id",
                   'logD', 'date_str', 'env_label',
                   "SF", "SB", "FS", "FB", "BS", "BF")]

# Keep only lines with phenotypes that are in the kinship matrix
pheno_subset <- pheno[pheno$pop_label %in% colnames(Ginv),]
pheno_subset <- pheno_subset[pheno_subset$env_label %in% c('NaCl', 'NGM'),]

# Scale covariates and traits
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")
pheno_subset$temperature <- as.numeric(pheno_subset$temperature)
pheno_subset$rel_humidity <- as.numeric(pheno_subset$rel_humidity)
pheno_subset$logD <- as.numeric(pheno_subset$logD)

for (j in c('temperature', "rel_humidity", "logD")) {
  pheno_subset[, j][is.na(pheno_subset[, j])] <- mean(pheno_subset[, j], na.rm = TRUE)
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE)) / sd(pheno_subset[, j], na.rm = TRUE)
}

for (j in vect_P_traits) {
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE)) / sd(pheno_subset[, j], na.rm = TRUE)
}

# Ensure the pop_label column is a character vector for proper merging
pheno_subset$pop_label <- as.character(pheno_subset$pop_label)
subset_genotype$pop_label <- as.character(subset_genotype$pop_label)

# Merge the data frames by pop_label
merged_data <- merge(pheno_subset, subset_genotype, by = 'pop_label')

# Plot phenotype data
pheno_long <- merged_data %>%
  pivot_longer(cols = c(FS, SF, BS, SB, BF, FB), names_to = "variable", values_to = "value")

pheno_long$variable <- factor(pheno_long$variable, levels = c("FS", "SF", "BS", "SB", "BF", "FB"))

# Calculate the mean of the traits
pheno_long_mean <- pheno_long %>%
  group_by(pop_label, variable) %>%
  summarize(value = mean(value, na.rm = TRUE))

# Initialize the final DataFrame for plotting
df_to_plot <- data.frame(pop_label = unique(pheno_long_mean$pop_label),
                         SF = rep(0, length(unique(pheno_long_mean$pop_label))),
                         SB = rep(0, length(unique(pheno_long_mean$pop_label))),
                         FS = rep(0, length(unique(pheno_long_mean$pop_label))),
                         FB = rep(0, length(unique(pheno_long_mean$pop_label))),
                         BS = rep(0, length(unique(pheno_long_mean$pop_label))),
                         BF = rep(0, length(unique(pheno_long_mean$pop_label))))

# Fill the DataFrame with the mean values
for(i in vect_P_traits){
  df_to_plot[, i] <- pheno_long_mean[pheno_long_mean$variable == i,]$value
}

# Merge the genotype data into the final DataFrame
df_to_plot <- merge(df_to_plot, merged_data[c('pop_label', 'Genotype')], by = 'pop_label')

plots <- list()

# Using the Set3 color palette from RColorBrewer
base_colors <- brewer.pal(9, "Set1")

# Augmenting the palette with additional colors
colors <- c(base_colors, rainbow(15 - length(base_colors)))

# Create scatter plots for each trait against SF
for(i in vect_P_traits){
  corr_coef <- cor(df_to_plot[["SF"]], df_to_plot[[i]], method = "pearson") # or method = "spearman"
  
  p <- ggplot(df_to_plot, aes_string(x = "SF", y = i, color = "Genotype")) +
    geom_point(show.legend = FALSE) +
    labs(x = 'SF', y = i) +
    ggtitle(paste("SF vs", i, "\nCorrelation:", round(corr_coef, 2))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  plots[[i]] <- p
}

# Arrange the plots in a grid layout and save as PDF
pdf(gsub('XXX', output, "Pheno_Trait_VS_traits_plot_XXX.pdf"), width = 9)
grid.arrange(grobs = plots, ncol = 5)
dev.off()
