### Script to plot SNPs along chromosomes find by GATK pipelines
## Input : *_snp_positions.css : can be snp_position for specific pop of all lines, just have to not have been hard filter
#          SNP_2exclude_All.tsv - file which contains SNPs name which hae to be remove
## Output : A plot wich summarize the lenght of each chomosome, the number of SNPs found on it, and the proportion corresponding

# Set the working directory
setwd("~/Documents/Worms/Plot_GATK/Pop")

# Load necessary libraries
library(ggplot2)
library(gridExtra)

# Load the data
All <- read.csv("A6_snp_positions.csv")
exclude <- read.csv('SNP_2exclude_All.tsv', header = FALSE, sep = '\t', col.names = colnames(All))

# Create a unique key for each SNP position
All$key <- paste(All$CHROM, All$POS, sep = '_')
exclude$key <- paste(exclude$CHROM, exclude$POS, sep = '_')

# Mark SNPs as 'FILTER' if they are in the exclude list, otherwise as 'PASS'
All$Filter <- as.factor(ifelse(All$key %in% exclude$key, 'FILTER', 'PASS'))

# Create a dummy variable for plotting
All$Plot <- rep(1, nrow(All))

# Filter the data to include only 'PASS' SNPs
All_filter <- All[All$Filter == 'PASS',]

# Plot the distribution of SNPs along chromosomes
ggplot(All_filter, aes(x = POS, y = Plot, color = Filter)) +
  geom_point(size = 0.5) +
  facet_wrap(~ CHROM, ncol = 2, scales = "free_x") +
  theme_bw() +
  ggtitle('SNP repartition for All populations along chromosomes')

# Calculate the number of SNPs and chromosome lengths
nbr <- c()
length_chrom <- c()
for(chrom in unique(All_filter$CHROM)){
  nbr <- c(nbr, nrow(All_filter[All_filter$CHROM == chrom,]))
  length_chrom <- c(length_chrom, All_filter[All_filter$CHROM == chrom,][nrow(All_filter[All_filter$CHROM == chrom,]),]$POS)
}

# Create a summary dataframe
Summary <- data.frame(Chromosome = unique(All$CHROM),
                      Number_of_SNPs = nbr,
                      Length_chromosome = length_chrom,
                      Proportion = round(nbr / length_chrom * 100, 5))

# Plot the length of chromosomes
g <- ggplot(Summary) +
  geom_point(aes(x = Chromosome, y = Length_chromosome, color = Chromosome)) +
  theme_minimal() +
  labs(y = 'Length of chromosome') +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.title.x = element_blank())

# Plot the number of SNPs per chromosome with proportions
g1 <- ggplot(Summary) +
  geom_bar(stat = "identity", aes(y = Number_of_SNPs, x = Chromosome, color = Chromosome, fill = Chromosome)) +
  geom_text(aes(y = Number_of_SNPs, x = Chromosome, label = paste(Proportion, '%')), vjust = -0.5, size = 3) +
  theme_minimal() +
  labs(y = 'Number of SNPs') +
  theme(legend.position = "none")

# Combine the plots
plots <- list(g, g1)

# Save the combined plots as a PDF
pdf('Distribution_SNP_pop.pdf')
grid.arrange(grobs = plots, ncol = 1, heights = c(1, 2))
dev.off()
