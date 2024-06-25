args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
# Name of output
output <- args[2]

# Load packages
library(vcfR)
library(ggplot2)
library(rutilstimflutre)

##### First step : kinship matrix ######

# Read vcf file
vcf <- read.vcfR(file)

# Extract genotype information
genotype_info <- vcf@gt

# Extract SNP positions
snp_positions <- vcf@fix[, c("CHROM", "POS")]



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

nbr_na <- data.frame(NA_values=rowSums(is.na(convert_genotype)))

g<-ggplot(nbr_na,aes(x = NA_values))+
  geom_histogram()+
  theme_bw()+
  xlim(c(1,20))+
  ggtitle('NA distribution per SNPs')

ggsave(filename = gsub('XXX',output,'XXX_NA_distribution.pdf'),plot=g)

snp_positions <- as.data.frame(snp_positions[complete.cases(convert_genotype),])
convert_genotype <- convert_genotype[complete.cases(convert_genotype),]


M <- t(convert_genotype)

# Function to clean the column names
clean_colnames <- function(colnames) {
  sapply(colnames, function(name) {
    parts <- strsplit(name, "_")[[1]]
    if (length(parts) > 1 && parts[1] == parts[2]) {
      return(parts[1])
    } else {
      return(name)
    }
  })
}

rownames(M) <- clean_colnames(rownames(M))
rownames(M) <- gsub('CeMee','',rownames(M))
rownames(M) <- gsub('_sorted','',rownames(M))
colnames(M) <- paste0(paste0(snp_positions$CHROM,'_'),snp_positions$POS)
## VanRaden 2008
VR <- estimGenRel(
  M,
  relationships = "additive",
  method = "vanraden1"
)

## noia Vitezica et al, 2017
noia <- estimGenRel(
  M,
  relationships = "additive",
  method = "noia"
)

calcular_similarite_genetique <- function(gmat) {
  m <- nrow(gmat) # markers
  gmat <- scale(t(gmat))
  gmat[is.na(gmat)] <- 0
  # u = PCs, v = marker loadings, d = singular values for each PC
  # so that xx = svdx$u %*% diag(svdx$d) %*% t(svdx$v)
  # from Patterson 2006, Hoffman 2013, K = XX'= US(US)'
  X <- tcrossprod(gmat) / m
  svdx <- svd(X)
  D <- diag(sqrt(svdx$d))
  US <- svdx$u %*% D
  K <- tcrossprod(US)
  K <- K / mean(diag(K))
  rownames(K) <- colnames(K) <- rownames(gmat)
  return(K)
}

A<-calcular_similarite_genetique(convert_genotype)

write.csv(snp_positions,gsub('XXX',output,'XXX_snp_positions.csv'),row.names = FALSE,quote = FALSE) #Need for backsolving
write.csv(M,gsub('XXX',output,'XXX_t_convert_genotype.csv'),quote = FALSE) #Need for backsolving
write.csv(VR,gsub('XXX',output,'Kinship_matrix_VanRaden_XXX.csv'),row.names = FALSE,quote = FALSE) 
write.csv(A,gsub('XXX',output,'Kinship_matrix_Luke_XXX.csv'),row.names = FALSE,quote = FALSE) 
write.csv(noia,gsub('XXX',output,'Kinship_matrix_noia_XXX.csv'),row.names = FALSE,quote = FALSE) 
