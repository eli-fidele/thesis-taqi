
#=================================================================================#
#                             RANDOM MATRIX FUNCTIONS
#=================================================================================#

RM_normal <- function(M, normal_args = c(0,1), symm = F){
  # Extract parameters
  mu <- normal_args[1]
  sd <- normal_args[2]
  # Create [M x M] matrix
  P <- matrix(rep(NA, M * M), ncol = M)  
  # Generate rows
  for(i in 1:M){P[i,] <- rnorm(n = M, mean = mu, sd = sd)}
  # Make symmetric if prompted
  if(symm == T){P <- equalize_triangles(P)}
  # Return the matrix
  P
}

# Generate random stochastic matrix of size M, with choice of row function {r_stochastic, r_zeros}
RM_stoch <- function(M, symm = F, sparsity = F, opt1 = T){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  if(sparsity){row_fn <- r_zeros} else {row_fn <- r_stochastic} # choose row function
  # Generate rows
  for(i in 1:M){P[i,] <- row_fn(M)}
  # Make symmetric (if prompted)
  if(symm == T){
    # Add Transpose
    if(opt1){P <- equalize_triangles(P)}
    # Retune the rows
    P <- retune_rows(P)
  }
  # Return the matrix
  P
}

#=================================================================================#
#                         SPECIAL RANDOM MATRIX FUNCTIONS
#=================================================================================#

# Generate a tridiagonal matrix with normal entries
RM_trid <- function(M, symm = F){
  diagonal <- rnorm(n = M, 0, 2)
  P <- diag(diagonal)
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- rnorm(n = M, 0, 1)
  # Return the matrix
  P
}

# (Erdos-Renyi Graph)
# p_sparse is a probability between [0,1) so edges are connected ~ Bern(p) 
RM_erdos <- function(M, p_sparse, stoch = F){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  for(i in 1:M){
    # generate current row
    curr_row <- runif(M,0,1)
    # sample number of zeros ~ Bin(n,o)
    num_zeros <- rbinom(1,M,p)
    choices <- sample(1:M, num_zeros) # Isomorphic to Erdos-Renyi graphs!
    curr_row[choices] <- 0
    if (stoch == T){
      curr_row <- curr_row/sum(curr_row)
      } # Normalize if to be stochastic
    # append the row
    P[i,] <- curr_row             
  }
  # Return the matrix
  P
}

#=================================================================================#
#                            STOCHASTIC ROW FUNCTIONS
#=================================================================================#

# generates stochastic rows of size M
r_stochastic <- function(M){
  prob <- runif(M,0,1)
  prob/sum(prob) # normalize
}

# generates same rows as in r_stochastic(M), but with introduced random sparsity
r_zeros <- function(M){
  prob <- runif(M,0,1)
  num_zeros <- sample(1:(M-1),1) # At most M-1 zeros, as to ensure stochastic property
  choices <- sample(1:M, num_zeros) # Choose edges to disconnect
  prob[choices] <- 0
  prob/sum(prob) # normalize
}

#=================================================================================#
#                         RANDOM MATRIX DIAGNOSTICS 
#=================================================================================#

# Obtain the nondiagonal entries of a row given its row index
nondiagonal_entries <- function(row, row_index){
  indices <- data.frame(idx = 1:length(row))
  indices <- indices %>% filter(idx != row_index)
  # return the row with the given indices
  row[as.numeric(indices[,])]
}

# Vectorize matrix entries to study their distribution
vectorize_matrix <- function(P){as.vector(P)}

# returns proportion of positive entries of any matrix P
pos_entries <- function(P){
  pos_entries <- length(matrix(P[P[,] > 0], nrow = 1))
  pos_entries/(length(P))   
}

#=========================================================================#
#                  (SYMMETRIC MATRICES) HELPER FUNCTIONS
#=========================================================================#

# Manually make equal the upper triangle and lower triangle of the matrix
equalize_triangles <- function(P){
  # Run over entry of the matrix
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      # Restrict view to one of the triangles (i < j): Lower Triangle
      if(i < j){P[i,j] <- P[j,i]} # Equalize lower and upper triangles
    }
  }
  P # Return Symmetric Matrix
}

#================================================================#
#              SYMMETRIC RANDOM MATRIX DIAGNOSTICS 
#================================================================#

# Check if a matrix is symmetric 
is_symmetric <- function(P){isSymmetric(P)}

normal_params <- function(entries){
  print(paste("Mean: ",round(mean(entries),3),sep=""))
  print(paste("Standard Deviation: ",round(sqrt(var(entries)),3),sep=""))
}

# For a given normal matrix, visualize its entries as a histogram
visualize_normal_entries <- function(P, normal_args){
  # Vectorize the matrix into a row vector of its entries
  elements_P <- data.frame(x = vectorize_matrix(P))
  # Extract parameters
  mu <- normal_args[1]
  sd <- normal_args[2]
  # Get theoretical distribution function
  normal_density <- function(x){dnorm(x, mean = mu, sd = sd)}
  # Plot the histogram of its entries
  entries_hist <- ggplot(data = elements_P, mapping = aes(x)) + 
    geom_histogram(bins = 20, aes(y = stat(density))) +
    stat_function(fun = normal_density)
  # Return plot
  entries_hist
}

#=========================================================================#
#                 STOCHASTIC RANDOM MATRIX DIAGNOSTICS 
#=========================================================================#

# Check if a matrix is stochastic
is_row_stochastic <- function(P){
  row_is_stoch <- rep(F, nrow(P))
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    row_is_stoch[i] <- (row_sum == 1)
  }
  !(F %in% row_is_stoch)
}

# Print out the sums of each of the rows (to check row-stochasticity)
row_sums <- function(P){
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    print(paste("Row ",i,": ",(row_sum),sep=""))
  }
}

# Print out the sums of each of the rows (to check row-stochasticity)
col_sums <- function(P){
  for(i in 1:ncol(P)){
    col_sum <- sum(P[,i])
    print(paste("Col ",i,": ",(col_sum),sep=""))
  }
}

#=========================================================================#
#                         HELPER FUNCTIONS 
#=========================================================================#

retune_rows <- function(P){
  M <- ncol(P)
  for(i in 1:M){
    # Start at M-1 index
    curr_idx <- M - i
    rows_counted <- (curr_idx + 1):M
    # Tune row w.r.t previously tuned rows
    P <- tune_row(P, curr_idx)
  }
  # Set diagonal such that rows sum to 1
  diag <- rep(0, ncol(P))
  for(i in 1:nrow(P)){
    row <- P[i, ]
    diag[i] <- (1 - sum(nondiagonal_entries(row, i)))
  }
  diag(P) <- diag
  P
}

tune_row <- function(P, curr_idx){
  cols_counted <- (1 + curr_idx):M
  #sum(P[i,rows_counted])
  for(i in 1:nrow(P)){
    #print(current_sum)
    allowed_sum <- 1 - sum(P[i, cols_counted])
    #print(allowed_sum)
    rem_cols <- 1:curr_idx
    P[i, rem_cols] <- P[i, rem_cols]*allowed_sum
    #print(P[i, rem_cols])
    P[rem_cols, i] <- t(P[i, rem_cols])
  }
  P
}
