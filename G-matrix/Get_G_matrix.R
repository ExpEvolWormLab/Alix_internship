### Script to get G-matrix from GBLUP

# Load necessary library
library(MCMCglmm)

# Define the name for the analysis
name <- 'VanRaden_A6_NaCl'

# Set the working directory
setwd("~/Documents/Worms/GBLUP/G-matrix")

# Generate the filename for the VCV file
VCV_file <- gsub('XXX', name, 'XXX_MCMCmodel_VCV.csv')

# Load the VCV file
VCV <- read.csv(VCV_file)

# Define the number of traits
nb_trait <- 6

# Compute the posterior mode of the VCV
post_dist = posterior.mode(VCV)

# Reshape the posterior distribution into a matrix
WI_G_mat = matrix(post_dist[1:nb_trait^2], nb_trait, nb_trait)

# Write the G-matrix to a CSV file
write.csv(WI_G_mat, gsub('XXX', name, 'XXX_G_matrix.csv'), quote = FALSE, row.names = FALSE)
