# Set working directory
setwd("~/Documents/Worms/VariantCalling/Hetero")

# Load necessary libraries
library(tidyr)
library(gridExtra)

# Load CSV files containing heterozygosity data
files <- list.files(pattern = "Heterozygosity_.*.csv")
populations <- sub("Heterozygosity_(.*).csv", '\\1', files)
dfs <- lapply(files, read.csv)

# Load HD regions from file
df <- read.delim("20220216_c_elegans_divergent_regions_strain.bed", header = FALSE)
df$V4 <- NULL  # Remove the fourth column
df <- unique(df)  # Keep only unique rows

# Function to merge overlapping intervals
merge_intervals <- function(df) {
  result <- list()  # Initialize an empty list to store results
  current_interval <- c(df$V2[1], df$V3[1])  # Start with the first interval
  
  for (i in 2:nrow(df)) {  # Loop through remaining intervals
    if (df$V2[i] <= current_interval[2]) {  # Check if intervals overlap
      current_interval[2] <- max(current_interval[2], df$V3[i])  # Merge intervals
    } else {
      result <- rbind(result, current_interval)  # Add non-overlapping interval to result
      current_interval <- c(df$V2[i], df$V3[i])  # Start new interval
    }
  }
  
  result <- rbind(result, current_interval)  # Add the last interval to result
  return(as.data.frame(result))  # Return result as data frame
}

# Store chrom5 results for each population
Df_chrom5 <- data.frame(X = integer(0),
                        Chromosome = character(0),
                        Bin_number = integer(0),
                        Start = integer(0),
                        End = integer(0),
                        Number_of_heterozygotes = numeric(0),
                        HD = numeric(0),
                        Bin_chrom = integer(0),
                        Population = character(0))

# Initialize list to store plots
plots <- list()

# Loop through each file and population
for (f in 1:length(files)) {
  pop <- populations[[f]]  # Get population name
  Results <- dfs[[f]]  # Get corresponding data frame
  Results[is.na(Results)] <- 0  # Replace NA values with 0
  
  # Apply function to merge intervals
  chrom <- unique(Results$Chromosome)  # Get unique chromosomes
  
  HD_df <- data.frame(Chromosome = character(0), Start = integer(0), End = integer(0))  # Initialize HD_df
  for (i in chrom) {  # Loop through each chromosome
    df_chrom <- df[df$V1 == i, ]  # Subset HD regions for current chromosome
    res <- merge_intervals(df_chrom)  # Merge overlapping intervals
    colnames(res) <- c("Start", "End")  # Rename columns
    res$Chromosome <- rep(i, nrow(res))  # Add chromosome name column
    HD_df <- rbind(HD_df, res)  # Combine results into HD_df
  }
  
  row.names(HD_df) <- NULL  # Remove row names
  
  # Initialize vectors to store results
  HD_r <- integer(length(Results$Chromosome))  # Vector for HD regions
  Bin_chrom <- integer(length(Results$Chromosome))  # Vector for bin numbers
  nbr <- 0  # Variable to keep track of the number of rows processed
  
  # Loop over each chromosome in the chrom vector
  for (c in chrom) {
    Results_chrom <- Results[Results$Chromosome == c, ]  # Subset Results for current chromosome
    HD_df_chrom <- HD_df[HD_df$Chromosome == c, ]  # Subset HD_df for current chromosome
    
    for (j in 1:nrow(HD_df_chrom)) {  # Loop over each row in HD_df_chrom
      start <- HD_df_chrom$Start[[j]]  # Get start position of HD region
      end <- HD_df_chrom$End[[j]]  # Get end position of HD region
      
      for (i in 1:nrow(Results_chrom)) {  # Loop over each row in Results_chrom
        Bin_chrom[i + nbr] <- i  # Store bin number
        
        # Check for overlap between HD region and Results bin
        if (start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
            start <= Results_chrom$Start[i] && end >= Results_chrom$End[i] || 
            start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
            start >= Results_chrom$Start[i] && end <= Results_chrom$End[i]) {
          HD_r[i + nbr] <- -0.05  # Mark overlap in HD_r
          res <- c(res, 1)  # Add 1 to res
        }
      }
    }
    nbr <- nbr + nrow(Results_chrom)  # Update number of rows processed
  }
  HD_r[HD_r == 0] <- NaN  # Set non-overlapping entries to NaN
  Results$HD <- HD_r  # Add HD column to Results
  Results$Bin_chrom <- Bin_chrom  # Add Bin_chrom column to Results
  
  # Convert Chromosome column to factor if it isn't already
  Results$Chromosome <- as.factor(Results$Chromosome)
  
  Chrom5 <- Results[Results$Chromosome == 'V', ]  # Subset Results for chromosome V
  Chrom5_zoom <- Chrom5  # Define zoomed region (currently same as Chrom5)
  
  Chrom5$Population <- rep(pop, nrow(Chrom5))  # Add Population column
  Df_chrom5 <- data.frame(rbind(Df_chrom5, Chrom5))  # Combine results into Df_chrom5
  
  # Create plot for current population
  g5 <- ggplot(Chrom5_zoom) +
    geom_point(aes(x = Start, y = Number_of_heterozygotes), alpha = 0.6, color = 'blue') +  # Plot heterozygosity
    geom_point(aes(x = Start, y = HD), color = "black", alpha = 0.6) +  # Plot HD regions
    ggtitle(pop) +  # Add title
    theme_bw() +  # Use theme_bw
    theme(axis.title.x = element_blank(),  # Remove x-axis title
          axis.text.x = element_blank(),  # Remove x-axis text
          axis.ticks.x = element_blank())  # Remove x-axis ticks
  
  plots[[f]] <- g5  # Store plot in list
}

# Create combined plot for all populations
ggplot(Df_chrom5) +
  geom_point(aes(x = Start, y = Number_of_heterozygotes, color = Population), alpha = 0.6) +  # Plot heterozygosity
  geom_point(aes(x = Start, y = HD), color = "black", alpha = 0.6) +  # Plot HD regions
  facet_wrap(~ Population, ncol = 3) +  # Facet by population
  theme_bw() +  # Use theme_bw
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        legend.position = "none")  # Remove legend

# Save combined plot to PDF
ggsave('Chromosome_5_Allpop.pdf')
