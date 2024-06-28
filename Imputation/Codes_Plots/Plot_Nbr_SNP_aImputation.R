### Script to compare the number of SNPs per population after imputation
## Input : .*.imputed_snp_positions.csv - generate with Get_position.R applied to imputed file when they are filter at this final step (hard-filtering)
## Output : Number_SNPs_aImputation.pdf - barplot of SNP number after imputation
#           Shared_SNPs_pop.pdf -  comparison of each pop to compute how many SNPs they shared, plot it with a barplot
#           Heatmap_shared_SNPs.pdf - heatmap which summarize all previous results

# Set the working directory
setwd("~/Documents/Worms/Plot_GATK/Pop")

# Load necessary libraries
library(tidyr)
library(ggplot2)
library(combinat)
library(dplyr)

# Load the CSV files
files <- list.files(pattern = ".*.imputed_snp_positions.csv")
files_name <- gsub('founders','F',sub("(.*).imputed_snp_positions.csv", '\\1', files))
dfs <- lapply(files, read.csv)

# Calculate the number of SNPs for each population
nbr <- c()
for(i in 1:length(files)){
  nbr <- c(nbr, nrow(dfs[[i]]))
}

# Create a dataframe for plotting
To_plot <- data.frame(nbr = nbr, pop = files_name)

# Plot the number of SNPs per population
ggplot(To_plot) +
  geom_bar(stat = "identity", aes(y = nbr, x = pop, color = pop, fill = pop)) +
  geom_text(aes(y = nbr, x = pop, label = nbr), vjust = -0.5) +
  theme_minimal() +
  labs(y = 'Number of SNPs') +
  theme(legend.position = "none")

# Save the plot
ggsave('Number_SNPs_aImputation.pdf')

# Generate pairwise combinations of populations
comb <- combn(1:length(files), 2)
comb_name <- combn(files_name, 2)

# Calculate the number of common SNPs between pairs of populations
nbr_common <- c()
name <- c()

for(i in 1:ncol(comb)){
  df1 <- dfs[[comb[, i][1]]]
  df2 <- dfs[[comb[, i][2]]]
  name <- c(name, paste(comb_name[, i][1], comb_name[, i][2], sep = ','))
  if(nrow(df1) == 0 || nrow(df2) == 0){
    nbr_common <- c(nbr_common, 0)
  } else {
    nbr_common <- c(nbr_common, nrow(intersect(df1, df2)))
  }
}

# Create a dataframe for the number of common SNPs
data <- data.frame(Comparison = name, Nbr_SNPs = nbr_common, nbr_string = sprintf("%.2e", round(nbr_common)))

# Remove zero values
data <- data[data$Nbr_SNPs != 0,]

# Plot the number of common SNPs between populations
ggplot(data, aes(x = Comparison, y = Nbr_SNPs, fill = Comparison)) +
  geom_bar(stat = "identity") +
  geom_point() +
  labs(x = "Traits", y = "Number of SNPs", title = "Number of SNPs Shared between Pop (EEV removed)") +
  theme_minimal() +
  geom_text(aes(x = Comparison, y = Nbr_SNPs, label = Nbr_SNPs), vjust = -0.5, size = 1.75) +
  scale_x_discrete(labels = function(x) gsub(",", "\n", x)) +
  theme(legend.position = "none", axis.text.x = element_text(size = 7))

# Save the plot
ggsave('Shared_SNPs_pop.pdf')

# Split population names for creating a heatmap matrix
data_split <- separate(data, Comparison, into = c("Pop1", "Pop2"), sep = ",")
data_split$Pop1 <- factor(data_split$Pop1, levels = c('F', 'A6', 'CA', 'GA', 'GM', 'GT', 'LR', 'SMR'))
data_split$Pop2 <- factor(data_split$Pop2, levels = c('F', 'A6', 'CA', 'GA', 'GM', 'GT', 'LR', 'SMR'))
data_split$nbr_string <- NULL

# Add diagonal values
diagonal_data <- data.frame(Pop1 = files_name, Pop2 = files_name, Nbr_SNPs = nbr)
diagonal_data <- diagonal_data[diagonal_data$Pop1 != 'EEV', ]

heatmap_data <- rbind(data_split, diagonal_data)
heatmap_data[heatmap_data$Pop2 == 'F' & heatmap_data$Pop1 == 'A6', ][, c('Pop1', 'Pop2')] <- c('F', 'A6')
heatmap_data[heatmap_data$Pop2 == 'F' & heatmap_data$Pop1 == 'CA', ][, c('Pop1', 'Pop2')] <- c('F', 'CA')

# Create the heatmap
ggplot(heatmap_data, aes(x = Pop1, y = Pop2, fill = Nbr_SNPs)) +
  geom_tile() +
  scale_fill_gradient(low = 'yellow', high = "red") +
  geom_text(aes(label = Nbr_SNPs), color = "black", size = 2) +
  theme_minimal() +
  labs(x = "Population 1", y = "Population 2", fill = "Number of SNPs", title = "Heatmap of Shared SNPs Between Populations") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the heatmap
ggsave('Heatmap_shared_SNPs.pdf')
