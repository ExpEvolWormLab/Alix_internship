### Script to plot univariate vs. multivariate SNP effects
## Outputs : Density_plots_comparions_Mt_Uni*pdf - compare density function of median of SNPs effect found in multivariate analysis and univariate analysis
#            Uni_VS_Multi*pdf - plot SNPs effect of MT analysis againt SNPs effect of ST analysis

# Set the working directory and define output directories
output <- 'VanRaden_A6_NaCl_0.99'
output_dir <- "~/Documents/Worms/GBLUP"
res_MT_dir <- "~/Documents/Worms/GBLUP/Pipeline_GBLUP/Pruned"
res_ST_dir <- "~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned"

setwd(res_MT_dir)
# Load the multivariate results CSV files
files <- list.files(pattern = gsub('XXX', output, "Summary_[A-Z]*_XXX.csv"))
motif_MT <- paste0("Summary_(.*)_", output, ".csv")
file_outputs <- sub(motif_MT, "\\1", files)
dfs_MT <- lapply(files, read.csv)

setwd(res_ST_dir)
# Load the univariate results CSV files
files <- list.files(pattern = gsub('XXX', output, "Summary_[A-Z]*_Uni_XXX.csv"))
motif_ST <- paste0("Summary_(.*)_Uni_", output, ".csv")
file_outputs <- sub(motif_ST, "\\1", files)
dfs_ST <- lapply(files, read.csv)

setwd(output_dir)

library(ggplot2)

# Combine the MT and ST data frames
combined_dfs <- data.frame()
combined_bis <- data.frame()
for (i in 1:length(dfs_MT)) {
  df_MT <- dfs_MT[[i]]
  df_ST <- dfs_ST[[i]]
  
  # Calculate confidence interval (IC)
  df_MT$IC <- abs(df_MT$lower - df_MT$upper)
  df_ST$IC <- abs(df_ST$lower - df_ST$upper)
  
  # Add source and trait information
  df_MT$source <- "MT"
  df_ST$source <- "ST"
  df_MT$Trait <- file_outputs[[i]]
  df_ST$Trait <- file_outputs[[i]]
  
  # Combine the data frames
  combined_df <- rbind(df_MT, df_ST)
  combined_dfs <- rbind(combined_df, combined_dfs)
  
  combined_bis <- rbind(combined_bis, merge(df_MT, df_ST, by = c('X', 'Trait', 'CHROM')))
}

# Rename columns for clarity
names(combined_bis)[names(combined_bis) == "median.x"] <- "median_MT"
names(combined_bis)[names(combined_bis) == "median.y"] <- "median_ST"

# Define colors for the plot
colors <- c("ST All" = "blue", "ST Credible" = "cyan", "MT All" = "pink", "MT Credible" = "red")

# Create the density plot with facets for each trait
ggplot() +
  geom_density(data = combined_dfs[combined_dfs$source == 'MT',], aes(x = median, fill = 'MT All'), color = 'black', alpha = 0.4) +
  geom_density(data = combined_dfs[combined_dfs$source == 'MT' & combined_dfs$credible == 'Credible',], aes(x = median, fill = 'MT Credible'), color = 'black', alpha = 0.4) +
  geom_density(data = combined_dfs[combined_dfs$source == 'ST',], aes(x = median, fill = 'ST All'), color = 'black', alpha = 0.4) +
  geom_density(data = combined_dfs[combined_dfs$source == 'ST' & combined_dfs$credible == 'Credible',], aes(x = median, fill = 'ST Credible'), color = 'black', alpha = 0.4) +
  labs(title = "Density of Median SNP Effects",
       x = "Median",
       y = "Density",
       fill = "Legend") +
  facet_wrap(~Trait, ncol = 1) +
  scale_fill_manual(values = colors) +
  expand_limits(x = range(combined_dfs$median)) +
  theme_minimal()

# Save the density plot
ggsave(gsub('XXX', output, 'Density_plots_comparions_Mt_Uni_XXX.pdf'))

# Add a new column for colors based on credibility
combined_bis$color <- with(combined_bis, ifelse(credible.y == 'Credible' & credible.x == 'Credible', 'Violet',
                                                ifelse(credible.y == 'Credible', 'Blue',
                                                       ifelse(credible.x == 'Credible', 'Red', 'Grey'))))

# Create the scatter plot with points colored by credibility
g <- ggplot(combined_bis, aes(x = median_ST, y = median_MT, color = color)) +
  geom_point(alpha = 0.6) +
  labs(title = "Plot of Median SNP Effects",
       x = "Median ST",
       y = "Median MT",
       color = "Credibility") +
  facet_wrap(~Trait, ncol = 2) +
  theme_minimal() +
  scale_color_manual(values = c('Blue' = 'blue', 'Red' = 'red', 'Violet' = 'purple', 'Grey' = 'grey'),
                     labels = c('Blue' = 'Credible in ST', 'Red' = 'Credible in MT', 'Violet' = 'Credible in Both', 'Grey' = 'Not Credible')) +
  theme(legend.position = 'right') +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black")

# Save the scatter plot
ggsave(gsub('XXX', output, 'Uni_VS_Multi_XXX.pdf'))
