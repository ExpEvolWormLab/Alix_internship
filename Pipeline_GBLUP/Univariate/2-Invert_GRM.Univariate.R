# Script to inverse the kinship matrix, they have to be inversible and the inverse has to be positive definite

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]
output <- args[2]

library(MASS)
library(Matrix)

A <- as.matrix(read.csv(file))

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

# Clean the column and row names
colnames(A) <- clean_colnames(colnames(A))
rownames(A) <- clean_colnames(rownames(A))

# Check the condition number
cond_number <- kappa(A)
tol <- 1e-10

if (cond_number > 1 / tol) {
  cat("Matrix is not invertible or is ill-conditioned (condition number:", cond_number, "). Computing pseudoinverse.\n")
  Ginv <- ginv(A)
} else {
  cat("Matrix is invertible (condition number:", cond_number, "). Computing inverse.\n")
  Ginv <- solve(A)
}

# Check if Ginv is positive definite
is.positive.definite <- function(x) {
  eigenvalues <- eigen(x, symmetric = TRUE, only.values = TRUE)$values
  all(eigenvalues > 0)
}

if (!is.positive.definite(Ginv)) {
  cat("Inverse is not positive definite. Adjusting to make it positive definite.\n")
  epsilon <- 1e-10
  diag(Ginv) <- diag(Ginv) + epsilon
  while (!is.positive.definite(Ginv)) {
    epsilon <- epsilon * 10
    diag(Ginv) <- diag(Ginv) + epsilon
  }
  cat("Adjusted epsilon to:", epsilon, "\n")
}

rownames(Ginv) <- colnames(A)
colnames(Ginv) <- colnames(A)
output_filename <- gsub(pattern = 'XXX', replacement = output, x = 'Inverted_kinship_matrix_XXX.csv')
write.csv(Ginv, output_filename, row.names = FALSE, quote = FALSE)
