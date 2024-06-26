### Script which modify phenotype table to be along with isotypes groups

## Load correpsondace table
correspondance_isotypes <- read.csv("~/Documents/Worms/VariantCalling/Isotypes/correspondance_isotypes.csv")

## Load phenotype table
pheno_table <- read.table("~/Documents/Worms/GBLUP/Transition_rates_estimates_may2024_export.txt",header = TRUE)

to_modify <- pheno_table$pop_label[which(pheno_table$pop_label %in% correspondance_isotypes$Line)]

pheno_table$pop_label[which(pheno_table$pop_label %in% correspondance_isotypes$Line)] <- correspondance_isotypes[correspondance_isotypes$Line %in% to_modify,]$Representant_line

write.csv(pheno_table,'~/Documents/Worms/GBLUP/Final_Transition_rates_estimates_may2024_export.csv',row.names = FALSE,quote = FALSE)
