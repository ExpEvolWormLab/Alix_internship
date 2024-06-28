#### Script to plot results of breeding value for univariate analysis
# Outputs : Trait_distribution_each_lines*pdf : breeding value distribution (median and sd) for each line along trait
#           forest_plot*pdf : depicts breeding values with are credible (interval of credibility of 95% doesn't overlap 0)
#           Trait_VS_traits_plot*pdf : plots breeding value for each traits against each other, compute Pearson coefficent to attest the correlation



# Setting the working directory
setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6")
output <- 'VanRaden_A6_NGM'
pheno_file <- 'Final_Transition_rates_estimates_may2024_export.csv'
M_file <- 'A6_t_convert_genotype.csv'


# Loading necessary libraries
library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)
library(HDInterval)
library(combinat)
library(RColorBrewer)

# Vector of traits
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")

# Reading in the genotype data
M <- as.data.frame(fread(M_file))
rownames(M) <- M[,1]  # Set row names
M[,1] <- NULL         # Remove the first column
M <- as.matrix(M)     # Convert to matrix

# Reading phenotype data
pheno <- read.csv(pheno_file)
# Keep only the lines which were phenotyped and genotyped
M <- M[rownames(M) %in% pheno$pop_label,]

### Trait breeding value mean for each line ###
# Initialize an empty dataframe for storing breeding values
hist_breeding_values <- data.frame(Line=character(0),
                                   Trait=character(0),
                                   Median=numeric(0),
                                   SD=numeric(0))


# Function to check if zero is within the interval
is_zero_in_interval <- function(lower, upper) {
  return(lower <= 0 & upper >= 0)
}

# Loop through each trait
motif <- paste0(output,"_TTT_MCMCmodel_Sol.csv")
for(p in vect_P_traits) {
  file <- gsub('TTT', p, motif)
  model_MCMC_WI_Sol <- read.csv(file)

  # Loop through each line
  for(name in rownames(M)) {
    # Select columns which contain the name of the line (i.e., breeding value distribution for each trait)
    col_subset <- colnames(model_MCMC_WI_Sol)[grepl(paste0(name, "$"), colnames(model_MCMC_WI_Sol))]
    if(length(col_subset) != 0) {
      subset <- data.frame(model_MCMC_WI_Sol[, col_subset]) # Consider just this line
      colnames(subset) <- gsub(paste0("^", "trait", "|", ".pop_label.", name, "$"), "", colnames(subset)) # Rename columns
      
      # Transform to a two-column dataframe
      sub_melt <- melt(subset)
      colnames(sub_melt) <- c('traits', 'breeding_value') # Rename columns
      
      
      # Calculate credible interval
      HDI <- hdi(subset) # 95%
      
      # Compute a summary for the line
      summary_df <- data.frame(
        Line = name,
        Trait = p,
        Median = colMedians(as.matrix(subset)),
        SD = colSds(as.matrix(subset)),
        Lower_ic = HDI["lower",],
        Upper_ic = HDI["upper",],
        Credible = !apply(HDI, 2, function(x) is_zero_in_interval(x["lower"], x["upper"])) # TRUE if 0 not in interval
      )
      
      # Add it to the existing data
      hist_breeding_values <- data.frame(rbind(hist_breeding_values, summary_df))
    }
  }
}

# Plot hist_breeding_values only for credible lines
Median_trait <- ggplot(hist_breeding_values, aes(x = Line, y = Median, color = Credible, group = Trait)) +
  geom_point(position = position_dodge(width = 0.9), size = 2) +
  geom_errorbar(aes(ymin = Median - SD, ymax = Median + SD), width = 0.2, linewidth = 0.25, position = position_dodge(width = 0.9), color = "black") +
  labs(x = "Lines", y = "Median", color = "Credible") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +  # Remove x-axis text
  facet_wrap(~Trait, ncol = 2, nrow = 3) +
  ggtitle(gsub('XXX', output, 'Breeding values XXX'))

Median_trait
pdf_filename <- gsub('XXX', output, "Trait_distribution_each_lines_XXX.pdf")
ggsave(pdf_filename, Median_trait)

##### Forest plot for credible lines
cred_BV <- hist_breeding_values[hist_breeding_values$Credible == TRUE,]
write.csv(cred_BV, gsub('XXX', output, 'Lines_credibles_XXX.csv'), row.names = FALSE, quote = FALSE)
length(unique(cred_BV$Line))

# Define colors for each Trait
colors <- c("SF" = "blue", "FS" = "green", "FB" = "red", "SB" = "cyan", "BS" = "pink", "BF" = "orange")

# Plot the forest plot
ggplot(cred_BV, aes(x = Median, y = Line, color = Trait)) +
  geom_point(position = position_dodge(width = 0.25), size = 2) +  # Points for the means with offset
  scale_color_manual(values = colors) +  # Apply the defined colors
  labs(title = "Forest Plot des Traits par Ligne",
       x = "Valeur Moyenne",
       y = "Ligne",
       color = "Trait") +
  theme_minimal() +  # Use a minimal theme
  theme(panel.grid.major = element_line(color = "gray"),
        panel.grid.minor = element_line(color = "gray"))

# Restore the default grid settings for future plots
theme_update(
  panel.grid.major = element_line(color = "gray90"),
  panel.grid.minor = element_line(color = "gray98")
)

ggsave(gsub('XXX', output, "forest_plot_XXX.pdf"), width = 4, height = 8)

### Plot to compare breeding value per trait: for each line, plot the value in one trait against another
df_to_plot <- data.frame(Line = unique(hist_breeding_values$Line),
                         SF = rep(0, nrow(hist_breeding_values)),
                         SB = rep(0, nrow(hist_breeding_values)),
                         FS = rep(0, nrow(hist_breeding_values)),
                         FB = rep(0, nrow(hist_breeding_values)),
                         BS = rep(0, nrow(hist_breeding_values)),
                         BF = rep(0, nrow(hist_breeding_values)))

to_compare <- unique(hist_breeding_values$Trait)
for(i in to_compare) {
  df_to_plot[, i] <- hist_breeding_values[hist_breeding_values$Trait == i,]$Median
}
comb <- combn(to_compare, 2)
plots <- list()

# Using the Set3 color palette from RColorBrewer
base_colors <- brewer.pal(9, "Set1")

# Augmenting the palette with additional colors
colors <- c(base_colors, rainbow(15 - length(base_colors)))

# Loop to create scatter plots comparing each pair of traits
for(i in 1:ncol(comb)) {
  x <- comb[, i][[1]]
  y <- comb[, i][[2]]
  
  corr_coef <- cor(df_to_plot[[x]], df_to_plot[[y]], method = "pearson") # Compute correlation
  
  p <- ggplot(df_to_plot, aes_string(x = x, y = y, color = factor(i))) +
    geom_point(show.legend = FALSE) +
    labs(x = x, y = y) +
    scale_color_manual(values = colors[i]) + 
    ggtitle(paste(x, "vs", y, "\nCorrelation:", round(corr_coef, 2))) +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank())
  
  plots[[i]] <- p
}

# Arrange all scatter plots in a grid layout
grid.arrange(grobs = plots, ncol = 5)
ggsave(gsub('XXX', output, "Trait_VS_traits_plot_XXX.pdf"))
