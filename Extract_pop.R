## Script which compute the distribution of population of the left 944
## Gentoyped left and genotyped and phenotyped left
setwd("~/Documents/Worms/Plot_GATK")

Df <- data.frame(Population = character(),
                 Num_generation = character(),
                 Num_line = character(),
                 Batch = character(),
                 stringsAsFactors = FALSE)

table <- read.csv('Final_lines.csv',header = FALSE)
colnames(table) <- c('Line')

### Function which extract pop, number of gen and line from the name of a line
## if a line is present in two batch it'll change the batch name by putting both of thme
extract_line <- function(table,Df,name_batch){
  Lines <- table[["Line"]]
  # Loop through the elements of Lines
  for (i in Lines) {
    if( !substr(i,1,2) %in% c("A6","CA","GA","GM","GT","LR","SM")){
      pop <- gsub("[0-9].*", "", i) #Get the first letter
      gen <- gsub("^[A-Za-z]+([0-9]+).*", "\\1",i) #Get number after the first letter
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch    
    }
    if (substr(i,1,2)=="A6"){
      pop <- substr(i, 1, 2)
      gen <- substr(i, 3, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch
    }
    if (substr(i,1,2) %in% c("CA","GA","GM","LR","GT")){
      pop <- substr(i, 1, 3)
      gen <- substr(i, 4, 5)
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch
    }
    if (substr(i,1,3)=="SMR"){
      pop <- substr(i, 1, 4)
      gen <- NaN
      n_line <- gsub(".*[A-Za-z](\\d+)$", "\\1", i) #Get number after the last letter
      batch <- name_batch 
    }
    if (substr(i,1,2)=="LR"){
      gen <- NaN
    }
    
    #print(paste("pop:", pop, "gen:", gen, "n_line:", n_line, "batch:", batch))
    # Create a new row with the extracted values
    matching_row <- Df$Population == pop & Df$Num_line == n_line & Df$Num_generation == gen
    #print(matching_row)
    if(any(matching_row)){
      Df$Batch[matching_row] <- paste0(Df$Batch[matching_row],paste0(", ",name_batch))
    } else {
      new_row <- data.frame(Population = pop,
                            Num_generation = gen,
                            Num_line = n_line,
                            Batch = batch)
      #print(new_row)
      # Add the new row to the existing data frame
      Df <- rbind(Df, new_row)
    }
  }
  return(Df)
}

Df <- extract_line(table,Df,"")
Df$Batch <- NULL

write.csv(table(Df$Population),'944_population_distribution.csv',row.names = FALSE,quote = FALSE)

pheno <- read.csv('Final_Transition_rates_estimates_may2024_export.csv')
table <- data.frame(Line=table[table$Line %in% pheno$pop_label,])

Df_pheno <- data.frame(Population = character(),
                 Num_generation = character(),
                 Num_line = character(),
                 Batch = character(),
                 stringsAsFactors = FALSE)
Df_pheno <- extract_line(table,Df_pheno,"")
Df_pheno$Batch <- NULL
write.csv(table(Df_pheno$Population),'944_population_distribution_phenotyped.csv',row.names = FALSE,quote = FALSE)
