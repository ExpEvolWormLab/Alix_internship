# Script to do the diagnostic of MCMC model for univariate analysis

# Capture command line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1] # Solution file (Sol)
file1 <- args[2] # Variance-covariance file (VCV)
session <- args[3] # Session identifier
trait <- args[4] # Trait to analyze
# Name of output file
output <- args[5]

# Load packages
options(rgl.useNULL=TRUE)
library(vcfR)
library(lattice)
library(MCMCglmm)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)

#### DIAGNOSTIC #### 

# Load the solution data
model_MCMC_WI_Sol <- read.csv(file)

# Load the variance-covariance data
model_MCMC_WI_VCV <- read.csv(file1)

# Define the variables to plot
to_plot <- c('logD', 'rel_humidity', 'temperature', trait)

# If session identifier is provided, include it in the variables to plot
if(session == 'session_id') {
  to_plot <- c('session_id', to_plot)
}

n <- length(to_plot) # Number of variables to plot

# Initialize a list to store the plots
plot_list <- list()

# Loop through each variable to plot autocorrelation
for(i in 1:n) { 
  x <- model_MCMC_WI_Sol[, i]  # Extract the column for the variable
  bacf <- acf(x, plot = FALSE)  # Compute autocorrelation function without plotting
  bacfdf <- with(bacf, data.frame(lag, acf))  # Convert ACF results to a data frame
  
  # Create the plot using ggplot2
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +  # Add horizontal line at y=0
    geom_segment(mapping = aes(xend = lag, yend = 0)) +  # Add segments for ACF values
    geom_hline(aes(yintercept = 0.05), linetype = "dotted", color = "red") +  # Add dotted lines for significance
    geom_hline(aes(yintercept = -0.05), linetype = "dotted", color = "red") +
    ylim(c(-0.5, 1)) +  # Set y-axis limits
    ggtitle(colnames(model_MCMC_WI_Sol)[i]) +  # Add title
    theme_bw() +  # Use theme_bw theme
    ylab('Correlation') +  # Label for y-axis
    xlab('Iteration') +  # Label for x-axis
    theme(plot.title = element_text(size = 11),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  
  # Add the plot to the list
  plot_list[[i]] <- q
}

# Loop through the first three columns of the VCV data for autocorrelation plots
for(j in 1:3) {
  x <- model_MCMC_WI_VCV[, j]  # Extract the column for the variable
  bacf <- acf(x, plot = FALSE)  # Compute autocorrelation function without plotting
  bacfdf <- with(bacf, data.frame(lag, acf))  # Convert ACF results to a data frame
  
  # Create the plot using ggplot2
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +  # Add horizontal line at y=0
    geom_segment(mapping = aes(xend = lag, yend = 0)) +  # Add segments for ACF values
    geom_hline(aes(yintercept = 0.05), linetype = "dotted", color = "red") +  # Add dotted lines for significance
    geom_hline(aes(yintercept = -0.05), linetype = "dotted", color = "red") +
    ylim(c(-0.5, 1)) +  # Set y-axis limits
    ggtitle(colnames(model_MCMC_WI_VCV)[j]) +  # Add title
    theme_bw() +  # Use theme_bw theme
    ylab('Correlation') +  # Label for y-axis
    xlab('Iteration') +  # Label for x-axis
    theme(plot.title = element_text(size = 11),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  
  # Add the plot to the list
  plot_list[[n + j]] <- q
}

# Save the plots to a PDF file
pdf(file = gsub("XXX", output, 'Model_pdf_MCMC_autocorr_XXX.pdf'))

# Number of plots per page
plot_per_page <- 12

# Loop through the plot list and save them in the PDF file
for (i in seq(1, length(plot_list), by = plot_per_page)) {
  end_index <- min(i + plot_per_page - 1, length(plot_list))  # Determine the end index for this page
  plot_list_subset <- plot_list[i:end_index]  # Extract the plots for this page
  
  # Print the plots in a grid with 3 columns
  print(grid.arrange(grobs = plot_list_subset, ncol = 3))
}

# Close the PDF device
dev.off()
