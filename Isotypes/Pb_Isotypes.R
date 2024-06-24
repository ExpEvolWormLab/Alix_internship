### Isotypes ###
# Script to create isotype group (Strain pairs with concordance > T)
# First and Second part should be split from others, and just produce the matrix in a .tsv file

# To change
# working directory
setwd("~/Documents/Worms/VariantCalling/Isotypes/")
# hard filtering file after imputation
file <- "I.final.hard_filter.vcf.gz"
# Matrix of concordance
matrix_name <- "matrix_concordance.by_hand.0.982.corrected.tsv"
nbr_na <- read.csv('Na_distribution.tsv')
# Cutoff
T <- 0.9945

library(ggplot2)
library(reshape2)

row.names(nbr_na) <- nbr_na$X
nbr_na$X <- NULL

# Write the matrix in .tsv file
matrix<-read.delim(matrix_name)

# Plot matrix distribution in order to set T
diag(matrix)<-NaN
df <- data.frame(melt(matrix))
ggplot(df, aes(x = value))+
  geom_histogram(fill = "skyblue", color = "black") +
  labs(title = "Concordance distribution", x = "Values", y = "Frequence")+
  geom_vline(xintercept = T,color='red')+
  theme_bw()

##### THIRD PART #####
### Select lines with a concordance higher than the T
### Create a list of potential friend : ie other line with it has more than T concordance
corr <- data.frame(matrix)
# Get lines which share more than T of similarity 
duo_isotype<-which(corr>T,arr.ind = TRUE)

# Extract row and column indices
row_indices <- duo_isotype[, 1]
col_indices <- duo_isotype[, 2]

# Create a data frame with row and column indices
scatter_data <- data.frame(x = row_indices, y = col_indices)

# Plot scatter plot
ggplot(scatter_data, aes(x = x, y = y)) +
  geom_point() +
  labs(title = "Scatter plot of duo_isotype indices", x = "Row indices", y = "Column indices")

name<-unique(sort(duo_isotype))
N_name <- length(name)

## Create a list of potential friend : for each line do list of other line with it has more than T concordanc
print('Find list of friends for each subset of lines which have concordance more than T')
start_time <- Sys.time()
liste_friend <- list()
n <- 1
for(k in name){ # Go through subset of lines with concordance > T
  progress_percentage <- (which(name==k) / length(name)) * 100
  cat(sprintf("\rProgress: %.2f%%", progress_percentage))
  print(colnames(matrix)[k])
  friends <- duo_isotype[,2][duo_isotype[,1]==k] # All lines it link to 
  fof <- c(k) # vector which contains lines link to k which are reciprocally link
  for(i in 1:(length(friends) - 1)) {
    progress_percentage1 <- (i / length(friends)) * 100
    cat(sprintf("\rProgress: %.2f%%", progress_percentage1))
    for(j in (i + 1):length(friends)) {
      if(all(c(friends[i], friends[j]) %in% duo_isotype)) { #if the two lines which has corr > T with k has corr > T between them
        fof <- c(friends[i], friends[j],fof) # we add them :)
      }
    }
  }
  liste_friend[[n]]<-unique(fof)
  n<-n+1
}
# Record the end time
end_time <- Sys.time()
# Calculate the time taken
time_taken <- end_time - start_time
cat(sprintf("End of list of friends in %.2f%%",time_taken))


##### FOURTH PART #####
### We're going to find subgroup in the previous group where all lines are linked with each other
### Function to do that in two part :
###     Part 1 : create a matrix which summarize relationship (friend or not friend)
###               find all the square around diag : obvious group of friends
###     Part 2 : compare the obvious group of friends with other in order to find hidden friend


find_true_friends <- function(liste,name,liste_friend) { 
  # Function to find group where they are all reciprocally linked (all friends !)
  # Put the relationship of our subset of line in a matrix
  # 1 is they have conco > T, 0 they haven't
  isotype_matrix <-
    matrix(0,
           nrow = length(liste),
           ncol = length(liste))
  # We look paiwise for concordance
  start_time <- Sys.time()
  print('Filling matrix')
  for (i in 1:(length(liste-1))) {
    num <- which(name == liste[i]) # find the index and so can have access to his friend list
    for (j in (i + 1):length(liste)) {
      if (liste[j] %in% liste_friend[[num]]) {
        # They are linked !
        isotype_matrix[i, j] <- 1
        isotype_matrix[i, j] <- 1
      }
    }
  }
  diag(isotype_matrix) <- 1 #They are linked to themselves
  
  # Record the end time
  end_time <- Sys.time()
  # Calculate the time taken
  time_taken <- end_time - start_time
  cat(sprintf("End of filling matrix in %.2f%%",time_taken))
  
  
  # Begin at the diagognal : want to find squares around diag
  isotype_groups <- list()
  g <- 1 #iterate list
  i <- 1 #line
  j <- 1 #col
  print('square finding')
  start_time <- Sys.time()
  while (i <= nrow(isotype_matrix)) { # While it's not the end of diagonal
    progress_percentage <- (i / nrow(isotype_matrix)) * 100
    cat(sprintf("\rProgress: %.2f%%", progress_percentage))
    flag <- c(1)
    while (length(unique(flag)) == 1 && # While every lines compared are linked (friend)
           j + k + 1 <= ncol(isotype_matrix)) { # While it's not the end of columns
      k <- k + 1 #number of test, test all cells at the right of diag until lines aren't linked
      for (n in 1:k) {
        if (isotype_matrix[i + n - 1, j + k] == 1) { # they are friends
          flag <- c(flag, 1)
        }
        else{ # they are not
          flag <- c(flag, 0)
        }
      }
      # Test are over
      if (length(unique(flag)) != 1) { # We're going to leave while loop, save the index of columns which were friends
        isotype_groups[[g]] <- i:(i + k - 1)
        g <- g + 1
      }
    }
    # Out of while loop
    if (length(unique(flag)) == 1) { # We leave the while loop because of the second condition
      #it's the end of the matrix
      isotype_groups[[g]] <- i:(nrow(isotype_matrix)) # Save which were friends
      i <- ncol(isotype_matrix) + 1
    }
    else{ # We leave the while loop because of the first condition
      if (k == 1) { # If it's failed at first test
        i <- i + 1
      }
      else{ # If it's failed after few tests
        i <- i + k - 1
      }
    }
    j <- i
  }
  # Record the end time
  end_time <- Sys.time()
  # Calculate the time taken
  time_taken <- end_time - start_time
  cat(sprintf("End of finding squares in %.2f%%",time_taken))
  
  # Now, we have list of groups to update because we could not have them all
  # Compare each groups to spot potential new friends' groups
  true_liste_isotypes_groups <- list() # list which store updated group
  united <- list() # list of group we unified
  gu <- 1
  g<-1
  print('Update groups')
  start_time <- Sys.time()
  if(length(isotype_groups)>1){ # Should not : to remove ?
    for(i in 2:(length(isotype_groups))){ # Go through previous group
      print(length(isotype_groups))
      for(j in (i-1):length(isotype_groups)-1){ # Compare pairwise
        group1<-isotype_groups[[i]]
        group2<-isotype_groups[[j]]
        
        flag<-1 # The friendship flag
        # Compare pairwise both group
        for(n1 in group1){
          for(n2 in group2){
            if(isotype_matrix[n1,n2]!=1){ # they are not friend !!
              flag<-0 # lower the friendship flag
              break # The ennemies of my friend are my ennemies : if one is not friend with the other, none are
            }
          }
          if(flag==0){ 
            break
          }
        }
        if(flag==1){ # Friends of my friends...
          ## To construct the true updated groups
          true_liste_isotypes_groups[[g]] <- union(group1,group2)
          g<-g+1
          united[[gu]] <- group1
          gu<-gu+1
          united[[gu]] <- group2
          gu<-gu+1
          # Update previous list in order to compare new groups with other
          isotype_groups[[length(isotype_groups)+1]]<- union(group1,group2)
        }
      }
    }
    # Updating group list
    to_add <- Filter(function(x) !any(sapply(united, function(y) any(x %in% y))), isotype_groups) # Add groups present in previous groups if they were not unified
    for(l in to_add){
      true_liste_isotypes_groups[[g]]<-l
      g<-g+1
    }
    
    filtered_list <- list()
  
    for (i in seq_along(true_liste_isotypes_groups)) {
      if (!any(sapply(true_liste_isotypes_groups[-i], function(x) all(true_liste_isotypes_groups[[i]] %in% x)))) { # removed group included in other
        filtered_list <- c(filtered_list, list(true_liste_isotypes_groups[[i]]))
      }
    }
    
    true_liste_isotypes_groups <- filtered_list
    
    # Record the end time
    end_time <- Sys.time()
    # Calculate the time taken
    time_taken <- end_time - start_time
    cat(sprintf("End of updating group in %.2f%%",time_taken))
    
    return(list(true_liste_isotypes_groups,isotype_matrix))
  }
  else{return(list(isotype_groups,isotype_matrix))}
}


liste_isotypes_groups <- c() # vector which is going to store the final isotype groups
liste_matrix <- list() # To follow what happens : to remove ?
i<-1
for (liste in liste_friend) {
  cat(sprintf("\rProgress: %.2f%%", i/length(liste_friend) * 100))
  liste_resultats <- find_true_friends(liste, name, liste_friend) # apply function
  
  isotypes_group <- lapply(liste_resultats[[1]], function(isotype_indices) { #convert name and order them
    sorted_isotypes <- sort(colnames(matrix)[liste[isotype_indices]])
    return(sorted_isotypes)
  })
  
  liste_isotypes_groups <- c(liste_isotypes_groups,isotypes_group)
  
  liste_matrix[[length(liste_matrix) + 1]] <- liste_resultats[[2]]
  
  i<-i+1
}

# Removed repeated and groups included in other
filtered_list <- list()

liste_isotypes_groups<-unique(liste_isotypes_groups)

for (i in seq_along(liste_isotypes_groups)) {
  if (!any(sapply(liste_isotypes_groups[-i], function(x) all(liste_isotypes_groups[[i]] %in% x)))) {
    filtered_list <- c(filtered_list, list(liste_isotypes_groups[[i]]))
  }
}

liste_isotypes_groups <- filtered_list

## Find isotypes which are shared between groups and remove them 
isotypes <- c()
for(i in liste_isotypes_groups){
  isotypes <- c(isotypes,i)
}
pb_isotype <- names(table(isotypes)[table(isotypes)!=1])

# Index des lignes à conserver
rows_to_keep <- rownames(matrix)[!(rownames(matrix) %in% pb_isotype)]

# Index des colonnes à conserver
cols_to_keep <- colnames(matrix)[!(colnames(matrix) %in% pb_isotype)]

new_matrix <- matrix[rows_to_keep, cols_to_keep]

write.table(new_matrix,'matrix_concordance.by_hand.0.982.corrected1.tsv',sep='\t',quote = FALSE)
write.csv(data.frame(Line=pb_isotype,Issue=rep('Shared_several_isotype_groups',length(pb_isotype))),'Removed_line_isotype_groups.csv',row.names = FALSE,quote = FALSE)

