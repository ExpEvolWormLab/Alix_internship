### Script to set cutoff for isotypes

setwd("~/Documents/Worms/VariantCalling/Isotypes/")

matrix <- read.csv('matrix_concordance.by_hand.0.982.tsv',sep = '\t')
doublons <- read.csv('Doublons.csv')

values <- c()
indice <- c()
indice1 <- c()
for(i in doublons$x){
  indice <- c(indice,i)
  indice1 <- c(indice1,paste0(i,'CeMee'))
  if(! i %in% rownames(matrix) || ! paste0(i,'CeMee') %in% rownames(matrix)){
    values <- c(values,NaN)
  }else{values <- c(values,matrix[i,paste0(i,'CeMee')])}
}

Df <- data.frame(i=indice,j=indice1,Concordance=values)

library(ggplot2)

# Assuming 'values' is your vector of data
# Create a data frame with the values
Plot <- data.frame(values)

# Tracé de la boîte à moustaches
ggplot(Df, aes(x = Concordance)) +
  geom_histogram(fill='red',color='black',binwidth = 0.003)+
  geom_hline(yintercept = 1)+
  geom_vline(xintercept = 0.9945)+
  theme_bw()

## Lines to removed
to_removed <- na.omit(Df$i[Df$Concordance<0.9945])
to_removed <- na.omit(c(to_removed,Df$j[Df$Concordance<0.9945]))


# Index des lignes à conserver
rows_to_keep <- !(rownames(matrix) %in% to_removed)

# Index des colonnes à conserver
cols_to_keep <- !(colnames(matrix) %in% to_removed)

new_matrix <- matrix[rows_to_keep, cols_to_keep]

write.table(new_matrix,'matrix_concordance.by_hand.0.982.corrected.tsv',sep='\t',quote = FALSE)
## To keep 
length(values[values>=0.99])


