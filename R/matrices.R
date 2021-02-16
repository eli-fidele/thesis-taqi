
#=================================================================================#
#                             RANDOM MATRIX FUNCTIONS
#=================================================================================#

RM_normal <- function(n, mean = 0, sd = 1, symm = F){
  # Create [n x n] matrix
  P <- matrix(rep(NA, n * n), ncol = n)  
  # Generate rows
  for(i in 1:n){P[i,] <- rnorm(n = n, mean, sd)}
  # Make symmetric if prompted
  if(symm == T){P <- make_symmetric(P)}
  # Return the matrix
  P
}

# Generate a tridiagonal matrix with normal entries
RM_trid <- function(n, symm = F){
  diagonal <- rnorm(n = n, 0, 2)
  P <- diag(diagonal)
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- rnorm(n = n, 0, 1)
  # Return the matrix
  P
}

# Generate random stochastic matrix of size n, with choice of row function {r_stochastic, r_zeros}
RM_stoch <- function(n, symm = F, sparsity = F){
  P <- matrix(rep(NA, n * n), ncol = n)  # create [n x n] transition matrix
  if(sparsity){row_fn <- r_zeros} else {row_fn <- r_stochastic} # choose row function
  # Generate rows
  for(i in 1:n){P[i,] <- row_fn(n)}
  # Make symmetric (if prompted)
  if(symm == T){
    # Equalize triangles
    P <- make_symmetric(P)
    # Nullify diagonal
    diag(P) <- rep(0, n)
    # Normalize rows
    for(i in 1:nrow(P)){
      row <- P[i, ]
      P[i,] <- row/sum(row)
    }
    # Set diagonal such that rows sum to 1
    diag <- rep(0, ncol(P))
    for(i in 1:nrow(P)){
      row <- P[i, ]
      diag[i] <- (1 - sum(nondiagonal_entries(row, i)))
    }
    diag(P) <- diag
    }
  # Return the matrix
  P
}

# (Erdos-Renyi Graph): p_sparse is a probability between [0,1) so edges are connected ~ Bern(p) 
RM_erdos <- function(n, p_sparse, stoch = F){
  P <- matrix(rep(NA, n * n), ncol = n)  # create [n x n] transition matrix
  p <- p_sparse # rename variable
  for(i in 1:n){
    # generate current row
    curr_row <- runif(n,0,1)
    # sample number of zeros ~ Bin(n,o)
    num_zeros <- rbinom(1,n,p)
    choices <- sample(1:n, num_zeros) # Isomorphic to Erdos-Renyi graphs!
    curr_row[choices] <- 0
    # Normalize if to be stochastic
    if (stoch == T){
      if(sum(curr_row) == 0){curr_row <- curr_row} else{
      curr_row <- curr_row/sum(curr_row)
      } 
    }
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
r_stochastic <- function(n){
  prob <- runif(n,0,1)
  prob/sum(prob) # normalize
}

# generates same rows as in r_stochastic(n), but with introduced random sparsity
r_zeros <- function(n){
  prob <- runif(n,0,1)
  num_zeros <- sample(1:(n-1),1) # At most n-1 zeros, as to ensure stochastic property
  choices <- sample(1:n, num_zeros) # Choose edges to disconnect
  prob[choices] <- 0
  prob/sum(prob) # normalize
}

#=========================================================================#
#                             HELPER FUNCTIONS
#=========================================================================#

# Manually make equal the upper triangle and lower triangle of the matrix
make_symmetric <- function(P){
  # Run over entry of the matrix
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      # Restrict view to one of the triangles (i < j): Lower Triangle
      if(i < j){P[i,j] <- P[j,i]} # Equalize lower and upper triangles
    }
  }
  P # Return Symmetric Matrix
}

# Obtain the nondiagonal entries of a row given its row index
nondiagonal_entries <- function(row, row_index){
  indices <- data.frame(idx = 1:length(row))
  indices <- indices %>% filter(idx != row_index)
  # return the row with the given indices
  row[as.numeric(indices[,])]
}

# Check if a matrix is stochastic
is_row_stochastic <- function(P){
  row_is_stoch <- rep(F, nrow(P))
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    row_is_stoch[i] <- (row_sum == 1)
  }
  !(F %in% row_is_stoch)
}

# returns proportion of positive entries of any matrix P
pos_entries <- function(P){
  pos_entries <- length(matrix(P[P[,] > 0], nrow = 1))
  pos_entries/(length(P))   
}
