### Isotypes ###
# Script to create isotype group (Strain pairs with concordance > 0.99985)
# First and Second part should be split from others, and just produce the matrix in a .tsv file

# To change
# working directory
setwd("~/Documents/Worms/VariantCalling/Isotypes")
# hard filtering file after imputation
file <- "I.final.hard_filter.vcf.gz"
# Cutoff
T <- 0.99985

library(vcfR)
library(ggplot2)

##### FIRST PART #####
### Read vcf file
### Convert genotype in number (-1 for Na, 0 or 2)

start_time <- Sys.time()

# Read vcf file
vcf <- read.vcfR(file)

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]

# Extract genotype information
genotype_info <- vcf@gt

# remove FORMAT column if it exists
if ("FORMAT" %in% colnames(genotype_info)) {
  genotype_info <- subset(genotype_info, select = -c(FORMAT))
}

# Function to convert genotypes
get_genotype <- function(genotype) {
  ifelse(substr(genotype, 1, 3) %in% c("0/0", "0|0"), 0,
         ifelse(substr(genotype, 1, 3) %in% c("1/1", "1|1"), 2,
                ifelse(substr(genotype, 1, 3) %in% c("./.", ".|."), NaN,NaN)))
}

convert_genotype <- apply(genotype_info, 2, get_genotype)

CeMeeV2<-read.csv('CeMeeV2_strains.txt',header = FALSE)
CeMeeV2$V1 <- paste0(CeMeeV2$V1,'CeMee')

col_subset <- unique(unlist(strsplit(as.character(CeMeeV2$V1), " ")))





nbr_na <- read.csv('Na_distribution.tsv')
rownames(nbr_na)<-nbr_na$X
nbr_na$X <- NULL

nbr_na_CeMee <- data.frame(X=rownames(nbr_na)[rownames(nbr_na) %in% col_subset],NA_values=nbr_na[rownames(nbr_na) %in% col_subset,])
rownames(nbr_na_CeMee)<-nbr_na_CeMee$X
nbr_na_CeMee$X <- NULL



g<-ggplot(nbr_na,aes(x=NA_values))+ 
  geom_histogram(fill='cyan',color='black',binwidth = 0.007)+ # each bin should have a width of 0.001 unit
  geom_vline(xintercept = 0.982,color='red')+
  labs(title = "Missing values distribution", x = "Values", y = "Number")+
  theme_bw()

ggsave('histo_NA.pdf',g)

g<-ggplot(nbr_na_CeMee,aes(x=NA_values))+ 
  geom_histogram(fill='cyan',color='black',binwidth = 0.007)+ # each bin should have a width of 0.001 unit
  geom_vline(xintercept = 0.982,color='red')+
  labs(title = "Missing values distribution", x = "Values", y = "Number")+
  theme_bw()

ggsave('histo_NA_CeMee.pdf',g)
  
to_keep <- rownames(nbr_na[nbr_na$NA_values<0.982,, drop = FALSE]) #drop prevent to remove rownames

for(i in rownames(nbr_na[nbr_na$NA_values>=0.982,, drop = FALSE])){ if(i %in% paste0(doublons$x,'CeMee')){print(i)}}

write.csv(data.frame(Line=rownames(nbr_na[nbr_na$NA_values>=0.982,, drop = FALSE]),Issue=rep('NA',length(rownames(nbr_na[nbr_na$NA_values>=0.982,, drop = FALSE])))),
        'Removed_line_NA.csv',
        quote = FALSE,
        row.names = FALSE)
reduced_convert_genotype <-convert_genotype[,to_keep]

##### SECOND PART #####
### Compute matrix which store how many SNPs are similar between two lines
#start_time <- Sys.time()
N_snp <- nrow(reduced_convert_genotype)
N_ind <- ncol(reduced_convert_genotype)
# Initialize matrix with zeros
matrix <- matrix(0, nrow = ncol(reduced_convert_genotype), ncol = ncol(reduced_convert_genotype))
k<-1
# Compare
for(i in 1:(N_ind)){
  progress_percentage <- (i / N_ind) * 100
  cat(sprintf("\rProgress: %.2f%%", progress_percentage))
  for (j in (i+1):N_ind-1){
    #matrix[i, j] <- sum(reduced_convert_genotype[, i] == reduced_convert_genotype[, j])
    subset <- na.omit(reduced_convert_genotype[,c(i,j)])
    matrix[i, j] <- sum(subset[, 1] == subset[, 2]) / nrow(subset)
    matrix[j, i] <- matrix[i, j] 
  }
}

colnames(matrix)<-colnames(reduced_convert_genotype)
rownames(matrix)<-colnames(reduced_convert_genotype)

# Record the end time
end_time <- Sys.time()

# Calculate the time taken
time_taken <- end_time - start_time

# Print the time taken
print(time_taken)

# Write the matrix in .tsv file
write.table(matrix,"matrix_concordance.by_hand.0.982.tsv",quote = FALSE,sep = '\t')
