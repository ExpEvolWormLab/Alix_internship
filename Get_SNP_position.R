## Script to get the position of SNPs present in a set of vcf file with similar names

setwd("~/Documents/Worms/Plot_GATK/Pop")

motif <- '(.*)_final.vcf.gz'
# Charger les fichiers vcf
files <- list.files(pattern = motif )

library(vcfR)

dfs <- lapply(files, read.vcfR)
files_name <- sub(motif ,"\\1",files)



##### First step : kinship matrix ######

for(i in dfs){
  vcf<-dfs[[i]]
  
  # Extract genotype information
  genotype_info <- vcf@gt
  
  # Extract SNP positions
  snp_positions <- vcf@fix[, c("CHROM", "POS")]
  
  write.csv(snp_positions,paste0(files_name,'_snp_positions.csv'),row.names = FALSE,quote = FALSE)
}

