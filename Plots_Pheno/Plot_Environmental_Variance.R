## Script to plot the environmentale variance for a specific population

# Set the working directory
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP")

# Define file paths and output name
matrix <- 'Inverted_kinship_matrix_VanRaden_A6.csv'
output <- 'A6_NaCl_NGM'
pheno_file <- 'Final_Transition_rates_estimates_may2024_export.csv'

# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Read the kinship matrix and clean column and row names
Ginv <- as.matrix(read.csv(matrix))
colnames(Ginv) <- gsub('CeMee', '', colnames(Ginv))
colnames(Ginv) <- gsub('_sorted', '', colnames(Ginv))
rownames(Ginv) <- colnames(Ginv)

# Read the phenotype data and select relevant columns
pheno <- read.csv(pheno_file)
pheno <- pheno[, c('pop_label', 'temperature', 'rel_humidity', "session_id",
                   'logD', 'date_str', 'env_label',
                   "SF", "SB", "FS", "FB", "BS", "BF")]

# Keep only lines with phenotypes and are in Ginv
pheno_subset <- pheno[pheno$pop_label %in% colnames(Ginv),]
pheno_subset <- pheno_subset[pheno_subset$env_label %in% c('NaCl', 'NGM'),]

# Scale covariates and reduce traits
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")

# Convert covariates to numeric and scale them
pheno_subset$temperature <- as.numeric(pheno_subset$temperature)
pheno_subset$rel_humidity <- as.numeric(pheno_subset$rel_humidity)
pheno_subset$logD <- as.numeric(pheno_subset$logD)

# Replace NA values with the mean and scale the covariates
for (j in c('temperature', "rel_humidity", "logD")) {
  pheno_subset[, j][is.na(pheno_subset[, j])] <- mean(pheno_subset[, j], na.rm = TRUE)
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE)) / sd(pheno_subset[, j], na.rm = TRUE)
}

# Scale the traits
for (j in vect_P_traits) {
  pheno_subset[, j] <- (pheno_subset[, j] - mean(pheno_subset[, j], na.rm = TRUE)) / sd(pheno_subset[, j], na.rm = TRUE)
}

# Calculate variance for each measure grouped by pop_label and env_label
variance_data <- summarise(group_by(pheno_subset, pop_label, env_label),
                           FB = var(FB, na.rm = TRUE),
                           BF = var(BF, na.rm = TRUE),
                           FS = var(FS, na.rm = TRUE),
                           SF = var(SF, na.rm = TRUE),
                           SB = var(SB, na.rm = TRUE),
                           BS = var(BS, na.rm = TRUE))

# Convert data from wide format to long format for plotting
variance_long <- gather(variance_data, key = "measure", value = "variance", -pop_label, -env_label)

# Plot the distribution of variances
ggplot(variance_long, aes(x = measure, y = variance, fill = env_label)) +
  geom_boxplot(position = position_dodge(0.9)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Measure", y = "Variance") +
  scale_fill_manual(values = c("NGM" = "darkorange", "NaCl" = "darkgreen")) +
  theme_bw()

# Save the plot
ggsave(gsub('XXX', output, 'Environmental_variance_XXX.pdf'))
