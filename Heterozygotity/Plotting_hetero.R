# Set working directory
setwd("~/Documents/Worms/VariantCalling/Hetero")

# Define results file and population of interest
Results_file <-  'Heterozygosity_GA.csv'
HD_file <- "20220216_c_elegans_divergent_regions_strain.bed"
pop <- 'A6'

# Load necessary libraries
library(ggplot2)

Results <- read.csv(Results_file)


# Replace NA values in Results with 0
Results[is.na(Results)] <- 0

# Load HD regions from file
df <- read.delim(HD_file, header = FALSE)
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

# Apply function to merge intervals for each chromosome
chrom <- unique(Results$Chromosome)  # Get unique chromosome names

HD_df <- data.frame(Chromosome = character(0), Start = integer(0), End = integer(0))  # Initialize empty data frame
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
      Bin_chrom[i+nbr] <- i  # Store bin number
      
      # Check for overlap between HD region and Results bin
      if (start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
          start <= Results_chrom$Start[i] && end >= Results_chrom$End[i] || 
          start <= Results_chrom$End[i] && end >= Results_chrom$Start[i] || 
          start >= Results_chrom$Start[i] && end <= Results_chrom$End[i]) {
        HD_r[i + nbr] <- -0.005  # Mark overlap in HD_r
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

# Create the plot
g <- ggplot(data = Results) +
  geom_point(aes(x = Bin_chrom, y = Number_of_heterozygotes, color = Chromosome), alpha = 0.6) +  # Plot heterozygosity
  geom_point(aes(x = Bin_chrom, y = HD), color = "black", alpha = 0.6) +  # Plot HD regions
  ggtitle(paste0("Heterozygosity along chromosome for ", pop)) +  # Add title
  theme_bw() +  # Use theme_bw
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank()) +  # Remove x-axis ticks
  facet_wrap(~ Chromosome, ncol = 2, scales = "free_x")  # Facet by chromosome

print(g)  # Print the plot
pdf_filename <- paste0("Heterozygotity_HDregion_", gsub(" ", "_", pop), ".pdf")  # Define PDF filename
ggsave(pdf_filename, g)  # Save plot to PDF

# Zoom in on interesting regions for Chromosome V
Chrom5 <- Results[Results$Chromosome == 'V', ]  # Subset Results for chromosome V
Chrom5_zoom <- Chrom5[Chrom5$Start > 16813367, ]  # Filter for start positions greater than 16813367
Chrom5_zoom <- Chrom5_zoom[Chrom5_zoom$End < 17627977, ]  # Filter for end positions less than 17627977

# Create zoomed-in plot
g5 <- ggplot(Chrom5_zoom) +
  geom_point(aes(x = Start, y = Number_of_heterozygotes), alpha = 0.6, color = 'blue') +  # Plot heterozygosity
  geom_point(aes(x = Start, y = HD), color = "black", alpha = 0.6) +  # Plot HD regions
  ggtitle(paste0("Heterozygosity along chromosome for ", pop)) +  # Add title
  theme_bw() +  # Use theme_bw
  theme(axis.title.x = element_blank(),  # Remove x-axis title
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks.x = element_blank())  # Remove x-axis ticks

pdf_filename <- paste0("Heterozygotity_HDregion_Zoom_V_", gsub(" ", "_", pop), ".pdf")  # Define PDF filename for zoomed plot
ggsave(pdf_filename, g5)  # Save zoomed plot to PDF

# Replace NA values with 0
Results[is.na(Results)] <- 0

# Modify specific HD values
Results[, Results$HD == -0.05] <- 1

# Define the number of repetitions for randomization
num_repetitions <- 50000

# Function for randomization
randomisation <- function(df, nbr) {
  df_random <- as.data.frame(df$HD)  # Create random dataframe
  colnames(df_random) <- "HD"  # Set column name
  df_random$Number_of_heterozygotes <- sample(df$Number_of_heterozygotes)  # Shuffle heterozygosity
  df_sorted <- df_random[order(-df_random$Number_of_heterozygotes), ]  # Sort by heterozygosity
  rownames(df_sorted) <- c(1:length(df_sorted$Number_of_heterozygotes))  # Set row names
  df_sorted <- df_sorted[1:nbr, ]  # Select top nbr rows
  prop_hetero <- sum(-df_sorted$HD) / nbr  # Calculate proportion of heterozygotes
  return(prop_hetero)  # Return proportion
}

# Initialize a list to store histograms
histograms <- list()

# Process each chromosome separately
for (chromosome in c('I','II','III','IV','V','X')) {
  chromosome_data <- Results[Results$Chromosome == chromosome, ]  # Filter Results for current chromosome
  
  Results_sorted <- chromosome_data[order(-chromosome_data$Number_of_heterozygotes), ]  # Sort by heterozygosity
  rownames(Results_sorted) <- c(1:length(Results_sorted$Chromosome))  # Set row names
  
  nbr <- 5 * length(Results_sorted$Chromosome) / 100  # Calculate top 5%
  Results_sorted <- Results_sorted[1:nbr, ]  # Select top 5%
  prop_hetero <- sum(-Results_sorted$HD) / nbr  # Calculate proportion of heterozygotes
  
  distribution <- replicate(num_repetitions, randomisation(chromosome_data, nbr))  # Randomize and calculate distribution
  
  hist <- ggplot(NULL, aes(x = distribution)) +
    geom_density(color = "black", fill = "white", bw = 0.00015) +  # Plot density
    geom_vline(xintercept = prop_hetero, color = "red") +  # Add vertical line at proportion
    ggtitle(paste("Chromosome", chromosome)) +  # Add title
    theme_bw()  # Use theme_bw
  
  histograms[[chromosome]] <- hist  # Store histogram
}

# Save histograms to a PDF file
pdf(paste0("Randomisation_", pop, ".pdf"))  # Open PDF file
grid.arrange(grobs = histograms)  # Arrange and save histograms
dev.off()  # Close PDF file

# Zoom in on Chromosome V and process similarly
Chrom5_zoom[is.na(Chrom5_zoom)] <- 0  # Replace NA with 0
Chrom5_zoom[, Chrom5_zoom$HD == -0.05] <- 1  # Modify specific HD values

Chrom5_zoom_sorted <- Chrom5_zoom[order(-Chrom5_zoom$Number_of_heterozygotes), ]  # Sort by heterozygosity
rownames(Chrom5_zoom_sorted) <- c(1:length(Chrom5_zoom$Chromosome))  # Set row names

nbr <- 5 * length(Chrom5_zoom$Chromosome) / 100  # Calculate top 5%
Chrom5_zoom_sorted <- Chrom5_zoom_sorted[1:nbr, ]  # Select top 5%
prop_hetero <- sum(-Chrom5_zoom_sorted$HD) / nbr  # Calculate proportion of heterozygotes

distribution_Chrom5 <- replicate(num_repetitions, randomisation(Chrom5_zoom, nbr))  # Randomize and calculate distribution

hist_Chrom5 <- ggplot(NULL, aes(x = distribution_Chrom5)) +
  geom_density(color = "black", fill = "white") +  # Plot density
  geom_vline(xintercept = prop_hetero, color = "red") +  # Add vertical line at proportion
  theme_bw()  # Use theme_bw

print(hist_Chrom5)  # Print histogram
pdf_filename <- paste0("Randomisation_Chrom5_", gsub(" ", "_", pop), ".pdf")  # Define PDF filename
ggsave(pdf_filename, hist_Chrom5)  # Save histogram to PDF
