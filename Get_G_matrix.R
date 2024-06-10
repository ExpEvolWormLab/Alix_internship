### Script to get G-matrix from GBLUP 
library(MCMCglmm)

name <- 'VanRaden_A6_NaCl'

setwd("~/Documents/Worms/GBLUP/G-matrix")

VCV_file <- gsub('XXX',name,'XXX_MCMCmodel_VCV.csv')

VCV <- read.csv(VCV_file)
nb_trait <- 6
post_dist = posterior.mode(VCV)
WI_G_mat=matrix(post_dist[1:nb_trait^2],nb_trait, nb_trait)

write.csv(WI_G_mat,gsub('XXX',name,'XXX_G_matrix.csv'),quote = FALSE,row.names = FALSE)


