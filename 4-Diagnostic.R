# To change
# Set working directory 
#setwd("~/Documents/Worms/GBLUP")
# Name of vcf file to use
#file <- 'founders.imputed.SNP.filtred.final.vcf.gz'
args <- commandArgs(trailingOnly = TRUE)
file <- args[1] # Sol
file1 <- args[2] # VCV
session <- args[3]
# Name of output
output <- args[4]
populations <- c("A6")
#populations <- c("A6",  "CA1",  "CA2",  "CA3",  "EEV",  "GA1",  "GA2",  "GA4",  "GM1",  "GM3",  "GT1",  "GT2",  "LR1",  "LR3", "SMR2", "SMR4")


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
model_MCMC_WI_Sol <- read.csv(file)
model_MCMC_WI_VCV <- read.csv(file1)
# Correlation plots

to_plot <- c('logD',
             'rel_humidity',
             'temperature',
             "SF",
             "SB",
             "FS",
             "FB",
             "BS",
             "BF")

if(session=='session_id'){
  to_plot <- c('session_id',to_plot)
}


n <- length(to_plot)

# Initialisation d'une liste pour stocker les graphiques
plot_list <- list()

for(i in 1:n){ ### Fixed effects and traits
  x <- model_MCMC_WI_Sol[,i]
  bacf <- acf(x, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    geom_hline(aes(yintercept = 0.05), linetype = "dotted", color = "red") +
    geom_hline(aes(yintercept = -0.05), linetype = "dotted", color = "red") +
    ylim(c(-0.5, 1)) +
    ggtitle(colnames(model_MCMC_WI_Sol)[i])+
    theme_bw()+
    ylab('Correlation')+
    xlab('Iteration')+
    theme(plot.title = element_text(size = 11),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  
  # Ajout du graphique à la liste
  plot_list[[i]] <- q
}


for(j in 1:36){
  x <- model_MCMC_WI_VCV[,j]
  bacf <- acf(x, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  q <- ggplot(data = bacfdf, mapping = aes(x = lag, y = acf)) +
    geom_hline(aes(yintercept = 0)) +
    geom_segment(mapping = aes(xend = lag, yend = 0)) +
    geom_hline(aes(yintercept = 0.05), linetype = "dotted", color = "red") +
    geom_hline(aes(yintercept = -0.05), linetype = "dotted", color = "red") +
    ylim(c(-0.5, 1)) +
    ggtitle(gsub('.pop_label','',colnames(model_MCMC_WI_VCV)[j]))+
    theme_bw()+
    ylab('Correlation')+
    xlab('Iteration')+
    theme(plot.title = element_text(size = 11),
          axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
  
  # Ajout du graphique à la liste
  plot_list[[n+j]] <- q
}

pdf(file = gsub("XXX",output,'Model_pdf_MCMC_autocorr_XXX.pdf'))
plot_per_page <- 12

for (i in seq(1, length(plot_list), by = plot_per_page)) {
  end_index <- min(i + plot_per_page - 1, length(plot_list))  # Déterminer l'index de fin pour cette page
  plot_list_subset <- plot_list[i:end_index]  # Extraire les graphiques pour cette page
  
  # Conversion de la liste en matrice pour l'affichage en 3 colonnes
  print(grid.arrange(grobs = plot_list_subset, ncol = 3))
}

dev.off()