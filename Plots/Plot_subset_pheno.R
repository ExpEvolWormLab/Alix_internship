### Script to plot phenotype data and compare traits using GBLUP

# Set the working directory
setwd("~/Documents/Worms/GBLUP")
setwd("~/Documents/Worms/GBLUP/Pipeline_GBLUP")

# Define file paths and output names
matrix <- 'Inverted_kinship_matrix_VanRaden_A6.csv'
output <- 'A6_NaCl_NGM'
pheno_file <- 'Final_Transition_rates_estimates_may2024_export.csv'

# Load necessary libraries
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Read the kinship matrix
Ginv <- as.matrix(read.csv(matrix))

# Correct column and row names
colnames(Ginv) <- gsub('CeMee', '', colnames(Ginv))
colnames(Ginv) <- gsub('_sorted', '', colnames(Ginv))
rownames(Ginv) <- colnames(Ginv)

# Read phenotype data
pheno <- read.csv(pheno_file)
pheno <- pheno[, c('pop_label', 'temperature', 'rel_humidity', "session_id",
                   'logD', 'date_str', 'env_label',
                   "SF", "SB", "FS", "FB", "BS", "BF")]

# Keep only lines with phenotypes that are in the kinship matrix
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

# Plot phenotype data
pheno_long <- pheno_subset %>%
  pivot_longer(cols = c(FS, SF, BS, SB, BF, FB), names_to = "variable", values_to = "value")

pheno_long$variable <- factor(pheno_long$variable, levels = c("FS", "SF", "BS", "SB", "BF", "FB"))

# Create the density plot with ggplot2
g <- ggplot(pheno_long, aes(x = value, fill = env_label, color = env_label)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~variable, scales = "free", ncol = 3) +
  theme_bw() +
  labs(title = "Trait Distribution by Condition",
       x = "Value",
       y = "Density") +
  theme(legend.position = "bottom")

# Save the density plot
ggsave(gsub('XXX', output, 'XXX_plots_pheno_subset.pdf'))

# Calculate the mean of the traits
pheno_long_mean <- pheno_long %>%
  group_by(pop_label, variable) %>%
  summarize(value = mean(value, na.rm = TRUE))

# Initialize the final DataFrame for plotting
df_to_plot <- data.frame(Line = unique(pheno_long_mean$pop_label),
                         SF = rep(0, length(unique(pheno_long_mean$pop_label))),
                         SB = rep(0, length(unique(pheno_long_mean$pop_label))),
                         FS = rep(0, length(unique(pheno_long_mean$pop_label))),
                         FB = rep(0, length(unique(pheno_long_mean$pop_label))),
                         BS = rep(0, length(unique(pheno_long_mean$pop_label))),
                         BF = rep(0, length(unique(pheno_long_mean$pop_label))))

# Fill the DataFrame with the mean values
for(i in vect_P_traits){
  df_to_plot[,i] <- pheno_long_mean[pheno_long_mean$variable == i,]$value
}

# Generate pairwise combinations of traits
comb <- combn(vect_P_traits, 2)
plots <- list()

# Using the Set3 color palette from RColorBrewer
base_colors <- brewer.pal(9, "Set1")

# Augmenting the palette with additional colors
colors <- c(base_colors, rainbow(15 - length(base_colors)))

# Create scatter plots for each pairwise combination of traits
for(i in 1:ncol(comb)){
  x <- comb[, i][[1]]
  y <- comb[, i][[2]]
  
  corr_coef <- cor(df_to_plot[[x]], df_to_plot[[y]], method = "pearson")
  
  p <- ggplot(df_to_plot, aes_string(x = x, y = y, color = factor(i))) +
    geom_point(show.legend = FALSE) +
    labs(x = x, y = y) +
    scale_color_manual(values = colors[i]) + 
    ggtitle(paste(x, "vs", y), subtitle = paste("Correlation:", round(corr_coef, 2))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  plots[[i]] <- p
}

# Save the scatter plots as a PDF
pdf(gsub('XXX', output, "Pheno_Trait_VS_traits_plot_XXX.pdf"), width = 9)
grid.arrange(grobs = plots, ncol = 5)
dev.off()
