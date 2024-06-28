### Script to compare the number of missing values (NAs) before and after imputation
# Input : .vchk for each vcf before and after imputation for each population, get .vchk by running analyse_vcf.sh
#      bash analyse_vcf.sh file.vcf.gz output_*a\bImputation 
# A specific motif is search, output format specified for analyse_vcf.sh need to be output_*[a or b]Imputation  
# Output : Impact_Imputation.pdf

# Set the working directory
setwd("~/Documents/Worms/Plot_GATK/Pop")

# Load necessary libraries
library(tidyr)
library(ggplot2)

# Load the CSV files for before and after imputation
files_b <- list.files(pattern = "output_.*_bImputation*")
files_name_b <- sub("output_(.*)_bImputation_PSC.vchk",'\\1',files_b)
dfs_b <- lapply(files_b, read.delim)

files_a <- list.files(pattern = "output_.*_aImputation*")
files_name_a <- sub("output_(.*)_aImputation_PSC.vchk",'\\1',files_a)
dfs_a <- lapply(files_a, read.delim)

# Initialize an empty dataframe to store the number of NAs
Nbr_Na <- data.frame(Population = character(0),
                     NA_b = numeric(0),
                     NA_a = numeric(0))

# Loop through each file to process the data
for(i in 1:length(files_a)){
  # Set columns name for after imputation
  colnames(dfs_a[[i]]) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions", "nIndels","average_depth","nSingleton","nHapRef","nHapAlt","nMissing")
  # Keep only the columns we are interested in
  dfs_a[[i]] <- dfs_a[[i]][c("sample","nRefHom","nNonRefHom","nHets","nSingleton","nMissing")]
  
  # Set columns name for before imputation
  colnames(dfs_b[[i]]) <- c("PSC","id","sample","nRefHom","nNonRefHom","nHets","nTransitions","nTransversions", "nIndels","average_depth","nSingleton","nHapRef","nHapAlt","nMissing")
  # Keep only the columns we are interested in
  dfs_b[[i]] <- dfs_b[[i]][c("sample","nRefHom","nNonRefHom","nHets","nSingleton","nMissing")]
  
  # Append the mean number of missing values to the Nbr_Na dataframe
  Nbr_Na <- data.frame(rbind(Nbr_Na, data.frame(Population = files_name_a[[i]],
                                         NA_b = mean(dfs_b[[i]]$nMissing),
                                         NA_a = mean(dfs_a[[i]]$nMissing))))
}

# Calculate the difference between the number of missing values before and after imputation
Nbr_Na$Difference <- Nbr_Na$NA_b - Nbr_Na$NA_a

# Create a plot with the distributions and differences
ggplot(data = Nbr_Na, aes(x = Population)) +
  geom_point(aes(y = NA_b, color = 'Before Imputation')) +
  geom_point(aes(y = NA_a, color = 'After Imputation')) +
  geom_text(data = subset(Nbr_Na, !is.na(Difference)), aes(y = NA_a, label = paste('+', round(Difference))), vjust = -1.5) +
  labs(y = "log(NA Values)", x = "Population", color = "Imputation") +
  scale_color_manual(values = c('Before Imputation' = 'blue', 'After Imputation' = 'red')) +
  theme_minimal() +
  theme(legend.position = "top") +
  scale_y_log10()  # Use a logarithmic scale

# Save the plot as a PDF
ggsave('Impact_Imputation.pdf')

