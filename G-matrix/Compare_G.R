## Script to compare different G-matrix

# Set working directory
setwd("~/Documents/Worms/GBLUP/G-matrix")

# Define the name for the analysis
name <- 'A6_NaCl'

# Generate the pattern to search for G-matrix files
to_seek <- gsub('XXX', name, ".*_XXX_G_matrix.csv")
# List files matching the pattern
files <- list.files(pattern = to_seek)
# Extract the names of the files
files_name <- sub(paste0('(.*)_', name, "_G_matrix.csv"), "\\1", files)
# Read the G-matrix files into a list
G_kin_matrices <- lapply(files, read.csv)

# Read the phenotypic G-matrix
Pheno <- read.csv(gsub('XXX', name, "XXX_MCMCmodel_VCV_G_matrix.csv"))

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(combinat)
library(gridExtra)

# Initialize a list to store the plots
plots <- list()

# Loop through each G-matrix file and create comparison plots
for (i in 1:length(files)) {
  # Melt the G-matrix data frame
  df1_melt <- melt(G_kin_matrices[[i]], varnames = c("Row", "Col"), value.name = "Value1")
  
  # Melt the phenotypic G-matrix data frame
  df2_melt <- melt(Pheno, varnames = c("Row", "Col"), value.name = "Value2")
  
  # Combine the two melted data frames by row and column indices
  combined_df <- data.frame(cbind(df1_melt, df2_melt))
  
  # Plot the values of G-matrix against the phenotypic matrix
  g <- ggplot(combined_df, aes(x = Value1, y = Value2)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # Add diagonal line
    ggtitle(paste0(files_name[i], ' VS Phenotypic')) +
    xlab(files_name[i]) +
    ylab('Phenotypic') +
    theme_bw()
  
  # Add the plot to the list
  plots[[i]] <- g
}

# Save the plots to a PDF file
pdf(gsub('XXX', name, "Comparison_G_matrix_XXX.pdf"))
grid.arrange(grobs = plots, ncol = 2)
dev.off()
