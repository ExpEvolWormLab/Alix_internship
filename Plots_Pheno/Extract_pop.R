## Script to compute the distribution of populations
## For both genotyped and genotyped + phenotyped lines

# Set working directory
setwd("~/Documents/Worms/Plot_GATK")

# Initialize an empty data frame
Df <- data.frame(Population = character(),
                 Num_generation = character(),
                 Num_line = character(),
                 Batch = character(),
                 stringsAsFactors = FALSE)

# Read the table containing the line information
table <- read.csv('Final_lines.csv', header = FALSE)
colnames(table) <- c('Line')

# Function to extract population, generation number, and line number from the line name
## If a line is present in two batches, it'll change the batch name by putting both of them
extract_line <- function(table, Df, name_batch) {
  Lines <- table[["Line"]]
  
  # Loop through the elements of Lines
  for (i in Lines) {
    if (!substr(i, 1, 2) %in% c("A6", "CA", "GA", "GM", "GT", "LR", "SM")) {
      pop <- gsub("[0-9].*", "", i) # Get the first letter
      gen <- gsub("^[A-Za-z]+([0-9]+).*", "\\1", i) # Get number after the first letter
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) # Get number after the last letter
      batch <- name_batch
    }
    if (substr(i, 1, 2) == "A6") {
      pop <- substr(i, 1, 2)
      gen <- substr(i, 3, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) # Get number after the last letter
      batch <- name_batch
    }
    if (substr(i, 1, 2) %in% c("CA", "GA", "GM", "LR", "GT")) {
      pop <- substr(i, 1, 3)
      gen <- substr(i, 4, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) # Get number after the last letter
      batch <- name_batch
    }
    if (substr(i, 1, 3) == "SMR") {
      pop <- substr(i, 1, 4)
      gen <- NaN
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) # Get number after the last letter
      batch <- name_batch
    }
    if (substr(i, 1, 2) == "LR") {
      gen <- NaN
    }
    
    # Check if the line is already present in the data frame
    matching_row <- Df$Population == pop & Df$Num_line == n_line & Df$Num_generation == gen
    
    if (any(matching_row)) {
      # If the line is already present, update the batch information
      Df$Batch[matching_row] <- paste0(Df$Batch[matching_row], paste0(", ", name_batch))
    } else {
      # Create a new row with the extracted values
      new_row <- data.frame(Population = pop,
                            Num_generation = gen,
                            Num_line = n_line,
                            Batch = batch)
      # Add the new row to the existing data frame
      Df <- rbind(Df, new_row)
    }
  }
  return(Df)
}

# Extract lines without batch information
Df <- extract_line(table, Df, "")
Df$Batch <- NULL

# Write the population distribution to a CSV file
write.csv(table(Df$Population), '944_population_distribution.csv', row.names = FALSE, quote = FALSE)

# Read the phenotype data
pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')

# Filter the table to include only lines present in the phenotype data
table <- data.frame(Line = table[table$Line %in% pheno$pop_label, ])

# Initialize an empty data frame for phenotyped lines
Df_pheno <- data.frame(Population = character(),
                       Num_generation = character(),
                       Num_line = character(),
                       Batch = character(),
                       stringsAsFactors = FALSE)

# Extract lines for phenotyped data
Df_pheno <- extract_line(table, Df_pheno, "")
Df_pheno$Batch <- NULL

# Write the population distribution for phenotyped lines to a CSV file
write.csv(table(Df_pheno$Population), '944_population_distribution_phenotyped.csv', row.names = FALSE, quote = FALSE)
