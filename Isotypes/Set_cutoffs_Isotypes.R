### Script to set cutoff for isotypes

# Set the working directory
setwd("~/Documents/Worms/VariantCalling/Isotypes/")

# Read the matrix of concordance values from a .tsv file
matrix <- read.csv('matrix_concordance.by_hand.0.982.tsv', sep = '\t')

# Read the list of duplicate entries from a .csv file
doublons <- read.csv('Doublons.csv')

# Define the concordance cutoff value
T <- 0.9945

# Initialize empty vectors to store values and indices
values <- c()
indice <- c()
indice1 <- c()

# Loop through each duplicate entry
for (i in doublons$x) {
  # Append the current entry to the indices
  indice <- c(indice, i)
  indice1 <- c(indice1, paste0(i, 'CeMee'))
  
  # Check if the current entry or its 'CeMee' version is not in the matrix
  if (!i %in% rownames(matrix) || !paste0(i, 'CeMee') %in% rownames(matrix)) {
    values <- c(values, NaN)
  } else {
    values <- c(values, matrix[i, paste0(i, 'CeMee')])
  }
}

# Create a data frame from the indices and concordance values
Df <- data.frame(i = indice, j = indice1, Concordance = values)

# Load the ggplot2 library for plotting
library(ggplot2)

# Assuming 'values' is your vector of data
# Create a data frame with the values for plotting
Plot <- data.frame(values)

# Plot the histogram of concordance values
ggplot(Df, aes(x = Concordance)) +
  geom_histogram(fill = 'red', color = 'black', binwidth = 0.003) + # Histogram with red fill and black borders
  geom_hline(yintercept = 1) + # Horizontal line at y = 1
  geom_vline(xintercept = T) + # Vertical line at the cutoff value
  theme_bw() # Apply a clean theme

# Identify lines to be removed based on the concordance cutoff
to_removed <- na.omit(Df$i[Df$Concordance < T])
to_removed <- na.omit(c(to_removed, Df$j[Df$Concordance < T]))

# Create an index for rows to keep
rows_to_keep <- !(rownames(matrix) %in% to_removed)

# Create an index for columns to keep
cols_to_keep <- !(colnames(matrix) %in% to_removed)

# Create a new matrix with the filtered rows and columns
new_matrix <- matrix[rows_to_keep, cols_to_keep]

# Write the new matrix to a .tsv file
write.table(new_matrix, 'matrix_concordance.by_hand.0.982.corrected.tsv', sep = '\t', quote = FALSE)
