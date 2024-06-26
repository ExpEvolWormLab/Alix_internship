## SINGLE TRAIT
#### RESULTS BREEDING VALUES AMONG PRUNING

# Set the working directory to the location of the results
setwd("~/Documents/Worms/GBLUP/Univariate_Pipeline_GBLUP/Results/Pruned")

## File and Variable Initialization
# Define the pattern for the file search
name <- 'VanRaden_A6_(0.[0-9]*)_NGM'
name1 <- 'VanRaden_A6_NGM'
# Define the genotype and phenotype file names
M_file <- 'A6_t_convert_genotype.csv'
pheno_file <- 'Final_Transition_rates_estimates_may2024_export.csv'

## Libraries and Setup
# Load necessary libraries for data manipulation and visualization
library(data.table)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)
library(HDInterval)
library(combinat)

# Create a pattern to search for specific CSV files
to_seek <- gsub('XXX',name,"XXX_MCMCmodel_Sol.csv")
# List all files matching the pattern
files <- list.files(pattern = to_seek)
# Extract the names from the files
files_name <- sub(to_seek, "\\1", files)
# Load each matching file into a list of data frames
dfs <- lapply(files, read.csv)
# Add an additional data frame to the list
dfs[[length(files)+1]] <- as.data.frame(fread(gsub('XXX',name1,"XXX_MCMCmodel_Sol.csv")))
# Append 'N' to the list of file names
files_name <- c(files_name, 'N')

# Read the genotype data file into a data frame
M <- as.data.frame(fread(M_file))
# Set the row names to the first column
rownames(M) <- M[,1]
# Remove the first column
M[,1] <- NULL
# Convert the data frame to a matrix
M <- as.matrix(M)

# Read the phenotype data file
pheno <- read.csv(pheno_file)
# Filter the genotype matrix to keep only the lines that are both phenotyped and genotyped
M <- M[rownames(M) %in% pheno$pop_label,]

## Trait Breeding Value Mean for Each Line
# Initialize a data frame to store breeding value statistics
hist_breeding_values <- data.frame(Line=character(0),
                                   Pruned=character(0),
                                   Trait=character(0),
                                   Median=numeric(0),
                                   SD=numeric(0),
                                   Lower_ic = numeric(0),
                                   Upper_ic = numeric(0),
                                   Credible = character(0))
# Initialize a list to store distribution of breeding values
distrib <- list() 

# Function to check if zero is within the credible interval
is_zero_in_interval <- function(lower, upper) {
  return(lower <= 0 & upper >= 0)
}

# Define the vector of trait names
vect_P_traits <- c("SF", "SB", "FS", "FB", "BS", "BF")

# Loop through each file in the list
for(k in 1:length(files_name)){
  model_MCMC_WI_Sol <- dfs[[k]] # Load the data frame for the current file
  i <- 1
  # Loop through each line in the genotype matrix
  for(name in rownames(M)){
    # Select columns that contain the name of the line (breeding value distribution for each trait)
    col_subset <- colnames(model_MCMC_WI_Sol)[grepl(paste0(name, "$"), colnames(model_MCMC_WI_Sol))]
    if(length(col_subset) != 0){ # If the line is present
      subset <- data.frame(model_MCMC_WI_Sol[, col_subset]) # Consider just this line
      # Rename columns to trait names
      colnames(subset) <- gsub(paste0("^", "trait", "|", ".pop_label.", name, "$"), "", colnames(subset))
      
      # Transform data into a two-column data frame
      sub_melt <- melt(subset)
      colnames(sub_melt) <- c('traits', 'breeding_value') # Rename columns
      i <- i + 1
      
      ## Find line with a credible breeding value
      # Calculate the 95% credible interval
      HDI <- hdi(subset)
      
      # Compute a summary for the line
      summary_df <- data.frame(
        Line = rep(name, 6),
        Pruned = files_name[[k]],
        Trait = vect_P_traits,
        Median = colMedians(as.matrix(subset)),
        SD = colSds(as.matrix(subset)),
        Lower_ic = HDI["lower",],
        Upper_ic = HDI["upper",],
        Credible = !apply(HDI, 2, function(x) is_zero_in_interval(x["lower"], x["upper"])) # TRUE if 0 not in interval
      )
      
      # Add the summary to the overall results data frame
      hist_breeding_values <- data.frame(rbind(hist_breeding_values, summary_df))
    }
  }
}

## Visualization
# Create a bar plot of the median breeding values along with pruning for each trait
ggplot(hist_breeding_values, aes(x=Pruned, y=Median, fill=Trait)) +
  geom_bar(stat = "identity") +
  labs(title = "Trait median along pruning",
       x = "Pruning",
       y = "Median") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3") +
  facet_wrap(~Trait, ncol = 2) +
  theme(legend.position = "none")

# Save the plot as a PDF file
ggsave('Median_BV_among_pruning.pdf')
