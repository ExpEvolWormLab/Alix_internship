### Compare G-matrix correlation of pruned data

# Set working directory
setwd("~/Documents/Worms/GBLUP/G-matrix")

# Define name and output for the analysis
name <- 'A6_NGM'
output <- 'A6_pruned_NGM'

# Pattern to search for G-matrix files
to_seek <- "VanRaden_(.*)_NGM_G_matrix.csv"
# List files matching the pattern
files <- list.files(pattern = to_seek)
# Extract names from the files
files_name <- sub(to_seek, "\\1", files)
# Read the G-matrix files into a list
G_kin_matrices <- lapply(files, read.csv)

# Add phenotypic G-matrix to the list
G_kin_matrices[[length(files) + 1]] <- read.csv(gsub('XXX', name, "XXX_MCMCmodel_VCV_G_matrix.csv"))
# Add 'Phenotypic' to the names
files_name <- c(files_name, 'Phenotypic')

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(combinat)
library(gridExtra)

# Create all combinations of the G-matrix indices
comb <- combn(1:length(files_name), 2)
comb_name <- combn(files_name, 2)
comb_name <- gsub('Luke', 'Hoffman', comb_name)

# Initialize vector to store correlation coefficients
corr_coef <- c()

# Calculate the Pearson correlation coefficient for each pair of G-matrices
for (i in 1:ncol(comb)) {
  df1_melt <- melt(G_kin_matrices[[comb[, i][1]]], varnames = c("Row", "Col"), value.name = "Value1")
  df2_melt <- melt(G_kin_matrices[[comb[, i][2]]], varnames = c("Row", "Col"), value.name = "Value2")
  corr_coef <- c(corr_coef, cor(df1_melt$Value1, df2_melt$Value2, method = "pearson"))
}

# Create a data frame for the correlation coefficients
comb_name_matrix <- apply(comb_name, 2, paste, collapse = " - ")
correlation_results <- data.frame(Pair = comb_name_matrix, Correlation = corr_coef)

# Print the correlation results
print(correlation_results)

# Create a correlation matrix
correlation_matrix <- matrix(NA, nrow = length(files_name), ncol = length(files_name))
rownames(correlation_matrix) <- files_name
colnames(correlation_matrix) <- files_name

# Fill the correlation matrix with calculated correlations
for (i in 1:ncol(comb)) {
  correlation_matrix[comb[1, i], comb[2, i]] <- corr_coef[i]
  correlation_matrix[comb[2, i], comb[1, i]] <- corr_coef[i]
}

# Fill the diagonal with 1s (perfect correlation with themselves)
diag(correlation_matrix) <- 1

# Melt the correlation matrix for ggplot2
melted_correlation_matrix <- melt(correlation_matrix, na.rm = TRUE)

# Plot heatmap
heatmap_plot <- ggplot(data = melted_correlation_matrix, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "purple", 
                       midpoint = 0.925, limit = c(0.85, 1), space = "Lab", 
                       name = "Pearson\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   size = 6, hjust = 1)) +
  coord_fixed()

# Save the heatmap plot
ggsave('Comparison_G_matrix_pruning.pdf')
