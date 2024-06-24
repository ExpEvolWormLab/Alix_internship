### Script to set cutoff for isotypes

setwd("~/Documents/Worms/VariantCalling/Isotypes/")

matrix <- read.csv('matrix_concordance.by_hand.0.982.tsv',sep = '\t')
doublons <- read.csv('Doublons.csv')
T <- 0.9945

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
  geom_vline(xintercept = T)+
  theme_bw()

## Lines to removed
to_removed <- na.omit(Df$i[Df$Concordance<T])
to_removed <- na.omit(c(to_removed,Df$j[Df$Concordance<T]))

### Write a file which contains all lines which have been removed and why
write.csv(data.frame(Line=to_removed,Issue=rep('LowConcordance',length(to_removed))),'Removed_line_concordance.csv',row.names = FALSE,quote = FALSE)

write.csv(data.frame(Line='CA150L45_sorted',Issue='Doublon'),'Removed_line_doublon.csv',row.names = FALSE,quote = FALSE)

to_removed<-c(to_removed,'CA150L45_sorted')


# Index des lignes à conserver
rows_to_keep <- !(rownames(matrix) %in% to_removed)

# Index des colonnes à conserver
cols_to_keep <- !(colnames(matrix) %in% to_removed)

new_matrix <- matrix[rows_to_keep, cols_to_keep]

write.table(new_matrix,'matrix_concordance.by_hand.corrected.tsv',sep='\t',quote = FALSE)
