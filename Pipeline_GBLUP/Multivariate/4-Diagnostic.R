## Script to diagnose MCMCglmm model - Multivariate analysis

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]  # Path to the Sol file
file1 <- args[2]  # Path to the VCV file
session <- args[3]  # Session identifier
output <- args[4]  # Output directory

# Load necessary packages
options(rgl.useNULL=TRUE)
library(vcfR)
library(lattice)
library(MCMCglmm)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(matrixStats)

#### DIAGNOSTIC #### 
# Load model solutions and variance-covariance matrices
model_MCMC_WI_Sol <- read.csv(file)
model_MCMC_WI_VCV <- read.csv(file1)

# Variables to include in the diagnostic plots
to_plot <- c('logD',
             'rel_humidity',
             'temperature',
             "SF",
             "SB",
             "FS",
             "FB",
             "BS",
             "BF")

# Include session_id if specified
if(session == 'session_id'){
  to_plot <- c('session_id', to_plot)
}

n <- length(to_plot)  # Number of variables to plot

# Initialize a list to store the plots
plot_list <- list()

# Loop to generate autocorrelation plots for fixed effects and traits
for(i in 1:n) {
  x <- model_MCMC_WI_Sol[, i]
  bacf <- acf(x, plot = FALSE)  # Calculate autocorrelation
  bacfdf <- with(bacf, data.frame(lag, acf))  # Convert to data frame
  
  # Generate the plot
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    geom_hline(aes(yintercept = 0.05), linetype = "dotted", color = "red") +
    geom_hline(aes(yintercept = -0.05), linetype = "dotted", color = "red") +
    ylim(c(-0.5, 1)) +
    ggtitle(colnames(model_MCMC_WI_Sol)[i]) +
    theme_bw() +
    ylab('Correlation') +
    xlab('Iteration') +
    theme(plot.title = element_text(size = 11),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  
  plot_list[[i]] <- q  # Add the plot to the list
}

# Loop to generate autocorrelation plots for variance-covariance components
for(j in 1:36) {
  x <- model_MCMC_WI_VCV[, j]
  bacf <- acf(x, plot = FALSE)  # Calculate autocorrelation
  bacfdf <- with(bacf, data.frame(lag, acf))  # Convert to data frame
  
  # Generate the plot
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    geom_hline(aes(yintercept = 0.05), linetype = "dotted", color = "red") +
    geom_hline(aes(yintercept = -0.05), linetype = "dotted", color = "red") +
    ylim(c(-0.5, 1)) +
    ggtitle(gsub('.pop_label', '', colnames(model_MCMC_WI_VCV)[j])) +
    theme_bw() +
    ylab('Correlation') +
    xlab('Iteration') +
    theme(plot.title = element_text(size = 11),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  
  plot_list[[n + j]] <- q  # Add the plot to the list
}

# Save all plots to a PDF
pdf(file = gsub("XXX", output, 'Model_pdf_MCMC_autocorr_XXX.pdf'))
plot_per_page <- 12  # Number of plots per page

# Loop to arrange and save the plots in the PDF
for (i in seq(1, length(plot_list), by = plot_per_page)) {
  end_index <- min(i + plot_per_page - 1, length(plot_list))  # Determine the end index for this page
  plot_list_subset <- plot_list[i:end_index]  # Extract plots for this page
  
  # Arrange the plots in a grid and print
  print(grid.arrange(grobs = plot_list_subset, ncol = 3))
}

dev.off()  # Close the PDF device
