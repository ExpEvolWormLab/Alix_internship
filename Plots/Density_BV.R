## Script to plot density function of breeding value in both analysis (MT, UT) in both condition (NGM, NaCl)
# Output : DensityPlot_BV*pdf 

# Define output name and file paths
output <- 'VanRaden_A6_NGM'
path_dir_Sol_MT <- "~/Documents/Worms/GBLUP/Pipeline_GBLUP/Results"
Sol_file <- 'VanRaden_A6_NGM_MCMCmodel_Sol.csv'
path_dir_Sol_ST <- "~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/A6"
M_file <- 'A6_t_convert_genotype.csv'
pheno_file <- 'Final_Transition_rates_estimates_may2024_export.csv'
output_dir <- "~/Documents/Worms/GBLUP"

# Librairies
library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)
library(HDInterval)
library(combinat)
library(RColorBrewer)


# Initialize Sol list to store data frames
Sol <- list()

# Set working directory and read the MT Sol file
setwd(path_dir_Sol_MT)
Sol[[1]] <- read.csv(Sol_file)
Sol_name <- c('Multi')

# Read genotype data
M <- as.data.frame(fread(M_file))
rownames(M) <- M[, 1]
M[, 1] <- NULL
M <- as.matrix(M)

# Set working directory and read univariate Sol files for each trait
setwd(path_dir_Sol_ST)
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")
for (p in vect_P_traits) {
  Sol_name <- c(Sol_name, 'Uni')
  Sol[[length(Sol) + 1]] <- read.csv(paste0(output, '_', p, '_MCMCmodel_Sol.csv'))
}

# Read phenotype data
pheno <- read.csv(pheno_file)
M <- M[rownames(M) %in% pheno$pop_label, ] # Keep only lines that were phenotyped and genotyped

# Set working directory to the root project folder
setwd(output_dir)

# Initialize data frame to store breeding value summaries
hist_breeding_values <- data.frame(Line = character(0),
                                   Trait = character(0),
                                   Median = numeric(0),
                                   Var = numeric(0),
                                   Lower_ic = numeric(0),
                                   Upper_ic = numeric(0),
                                   Credible = character(0),
                                   source = character(0))

# Function to check if zero is within the credible interval
is_zero_in_interval <- function(lower, upper) {
  return(lower <= 0 & upper >= 0)
}

# Iterate through each model's Sol data frame
i <- 1
for (model_MCMC_WI_Sol in Sol) {
  for (name in rownames(M)) {
    # Select columns containing the name of the line (i.e., breeding value distribution for each trait)
    col_subset <- colnames(model_MCMC_WI_Sol)[grepl(paste0(name, "$"), colnames(model_MCMC_WI_Sol))]
    if (length(col_subset) != 0) { # If the line is present
      subset <- data.frame(model_MCMC_WI_Sol[, col_subset]) # Consider just this line
      colnames(subset) <- gsub(paste0("^", "trait", "|", ".pop_label.", name, "$"), "", colnames(subset)) # Rename with trait labels
      
      # Calculate credible interval
      HDI <- hdi(subset) # 95%
      
      if (length(sub('trait([A-Z]*).pop_label.*', '\\1', col_subset)) != 1) {
        trait <- sub('trait([A-Z]*).pop_label.*', '\\1', col_subset)
      } else {
        trait <- vect_P_traits[[i - 1]]
      }
      
      # Compute a summary for the line
      summary_df <- data.frame(
        Line = rep(name, length(col_subset)),
        Trait = trait,
        Median = colMedians(as.matrix(subset)),
        Sd = colSds(as.matrix(subset)),
        Lower_ic = HDI["lower", ],
        Upper_ic = HDI["upper", ],
        Credible = !apply(HDI, 2, function(x) is_zero_in_interval(x["lower"], x["upper"])), # TRUE if 0 not in interval
        source = rep(Sol_name[i], length(col_subset))
      )
      
      # Add summary to the main data frame
      hist_breeding_values <- data.frame(rbind(hist_breeding_values, summary_df))
    }
  }
  i <- i + 1
}

# Define colors for each category with new labels
colors <- c("MT All" = "blue", "MT Credible" = "cyan", "ST All" = "pink", "ST Credible" = "red")

# Create density plot for median breeding values with normalization and new labels
m <- ggplot() +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Multi', ], aes(x = Median, fill = 'MT All'), color = 'black', alpha = 0.4) +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Uni', ], aes(x = Median, fill = 'ST All'), color = 'black', alpha = 0.4) +
  labs(title = "Density of Median Breeding values",
       x = "Median",
       y = "Density",
       fill = "Legend") +
  scale_fill_manual(values = colors) +
  theme_minimal()

# Create density plot for standard deviations of breeding values with normalization and new labels
s <- ggplot() +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Multi', ], aes(x = Sd, fill = 'MT All'), color = 'black', alpha = 0.4) +
  geom_density(data = hist_breeding_values[hist_breeding_values$source == 'Uni', ], aes(x = Sd, fill = 'ST All'), color = 'black', alpha = 0.4) +
  labs(title = "Density of Sd breeding values",
       x = "Sd",
       y = "Density",
       fill = "Legend") +
  scale_fill_manual(values = colors) +
  theme_minimal()

# Save the combined density plot to a PDF file
ggsave(gsub('XXX', output, 'DensityPlot_BV_XXX.pdf'))

# Arrange and display the two density plots side by side
grid.arrange(m, s, nrow = 1)
