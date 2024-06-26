### Script which modify phenotype table to be along with isotypes groups

## To change
# Set working directory 
setwd("~/Documents/Worms/VariantCalling/Isotypes")
# Name of correspondance file
correspondance_file <- "correspondance_isotypes.csv"
# Name of phenotype table to modify (need to have "pop_label" column)
pheno_file <- "~/Documents/Worms/GBLUP/Transition_rates_estimates_may2024_export.txt"
# Name of output 
output <- "Final_Transition_rates_estimates_may2024_export.csv"

## Load correpsondace table
correspondance_isotypes <- read.csv(correspondance_file)

## Load phenotype table
pheno_table <- read.table(pheno_file,header = TRUE)

to_modify <- pheno_table$pop_label[which(pheno_table$pop_label %in% correspondance_isotypes$Line)]

pheno_table$pop_label[which(pheno_table$pop_label %in% correspondance_isotypes$Line)] <- correspondance_isotypes[correspondance_isotypes$Line %in% to_modify,]$Representant_line

write.csv(pheno_table,output,row.names = FALSE,quote = FALSE)
