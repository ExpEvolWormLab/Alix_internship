### Script to analyze SNP effects for univariate analysis
## Outputs: 
## 1. Manhattan Plot - median of posterior SNP effects along chromosome colored by credibility
## 2. Density Plot - density function of median of posterior distribution of SNP effects

# Define output identifier and set working directory
output <- 'Uni_VanRaden_A6_NaCl_0.99'
setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")

# Load necessary libraries
library(ggplot2)
library(matrixStats)
library(tidyr)
library(gridExtra)

# Load CSV files for summary statistics
# Use pattern matching to find relevant files and read them into a list of data frames
files <- list.files(pattern = gsub('XXX', output, "Summary_[A-Z]*_XXX.csv"))
file_outputs <- sub(paste0("Summary_(.*)_", output, ".csv"), "\\1", files)
dfs <- lapply(files, read.csv)

# Initialize a list to store plots
plots <- list()

# Load SNP count file and keep only SNPs associated with one trait
snp_df <- read.csv(gsub('XXX', output, "SNP_Count_XXX.csv"))
snp_df <- snp_df[snp_df$Count == 1,]

# Loop through each summary statistics file
for(i in 1:length(dfs)){
  # Extract and preprocess data
  df <- dfs[[i]]
  df$credible <- factor(df$credible, levels = c("Not Credible", "Credible"))  # Convert credible column to factor
  df <- separate(df, X, into = c("CHROM", "POS"), sep = "_", remove = FALSE)  # Split the SNP identifier into CHROM and POS
  df$POS <- as.numeric(df$POS)  # Convert POS to numeric
  
  # Get the trait name from the file output list
  trait <- file_outputs[[i]]
  
  # Filter for SNPs associated only with this trait
  df_1 <- df[df$X %in% snp_df$SNP,]
  
  # Create a Manhattan plot for the current trait
  manhattan_plot <- ggplot(df_1, aes(x = POS, y = median, color = credible, alpha = credible)) +
    geom_point(size = 0.5) + 
    scale_color_manual(values = c("Not Credible" = 'grey', "Credible" = "red")) +  # Set colors for credible and not credible SNPs
    scale_alpha_manual(values = c("Not Credible" = 0.01, "Credible" = 1)) +  # Set transparency for credible and not credible SNPs
    theme_minimal() +  # Use a minimal theme for the plot
    ggtitle(gsub('TTT', trait, "Manhattan Plot Bayésien TTT")) +
    labs(x = "SNPs",
         y = "Effet médian à posteriori") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +  # Remove x-axis text and ticks
    facet_wrap(~ CHROM, ncol = 2, scales = "free_x")  # Create facets for each chromosome
  
  # Save the Manhattan plot as a PDF file
  ggsave(gsub('TTT', trait, gsub('XXX', output, 'Manhattan_plot_XXX_TTT.pdf')))
  
  # Filter for credible SNPs for density plot
  df_cred <- df[df$credible == 'Credible',]
  
  # Create a density plot for the current trait
  density_plot <- ggplot() +
    geom_density(data = df, aes(x = median), color = 'blue', fill = 'blue', alpha = 0.3) +  # Density plot for all SNPs
    geom_density(data = df_cred, aes(x = median), color = 'red', fill = 'red', alpha = 0.3) +  # Density plot for credible SNPs
    labs(title = paste("Density median SNPs effects -", trait),
         x = "Médiane",
         y = "Densité") +
    theme_minimal()  # Use a minimal theme for the plot
  
  # Add the density plot to the list of plots
  plots[[i]] <- density_plot
}

# Save all density plots as a single PDF file
pdf(gsub('XXX', output, 'Density_plot_XXX.pdf'))
grid.arrange(grobs = plots)
dev.off()
