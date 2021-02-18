
#=================================================================================#
#                           NORMAL RANDOM MATRICES 
#=================================================================================#

RM_normal <- function(N, mean = 0, sd = 1, symm = F, complex = F, hermitian = F){
  # Create [n x n] matrix with normally distributed entries
  P <- matrix(rnorm(N^2, mean, sd), nrow = N)  
  # Make symmetric if prompted
  if(symm || hermitian){P <- make_hermitian(P)}
  # Returns a matrix with complex entries if prompted
  if(complex){
    if(hermitian){
      P <- P + 1i * RM_normal(N, mean, sd, symm = T)
    } else{
      P <- P + 1i * RM_normal(N, mean, sd, symm = F)
    }
  }
  P <- P/sqrt(N) # Rescale the matrix
  P # Return the matrix
}

# Generate a tridiagonal matrix with normal entries
RM_trid <- function(N, symm = F){
  diagonal <- rnorm(N = N, 0, 2)
  P <- diag(diagonal)
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- rnorm(n = n, 0, 1)
  P# Return the matrix
}

# Generate a Gaussian (Hermite) Beta Ensemble matrix with Non-Invariant Dumitriu's Tridiagonal Model
RM_beta <- function(N, beta, complex = F){
  # Set the diagonal as a N(0,2) distributed row.
  diagonal <- rnorm(N, mean = 0, sd = 2)
  P <- diag(diagonal)
  # Set the off-1 diagonals as chi squared variables with df(beta), as given in Dumitriu's model
  df_seq <- beta*(N - seq(1,N-1)) # Get degrees of freedom sequence for offdigonal
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- rchisq(N-1, df_seq) # Generate tridiagonal
  # Add complex entries, if prompted
  if(complex){P <- P + (1i * RM_beta(N, beta))}
  P <- P/sqrt(2) # Rescale the entries by 1/sqrt(2)
  P # Return the matrix
}

#=================================================================================#
#                           STOCHASTIC RANDOM MATRICES 
#=================================================================================#

# Generate random stochastic matrix of size n, with choice of row function {r_stochastic, r_zeros}
RM_stoch <- function(N, symm = F, sparsity = F){
  P <- matrix(rep(NA, N * N), ncol = N)  # create [N x N] transition matrix
  if(sparsity){row_fn <- r_zeros} else {row_fn <- r_stochastic} # choose row function
  # Generate rows
  for(i in 1:N){P[i,] <- row_fn(N)}
  # Make symmetric (if prompted)
  if(symm == T){
    # Make lower and upper triangles equal
    P <- make_hermitian(P)
    # Nullify diagonal
    diag(P) <- rep(0, N)
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

# An Erdos-Renyi Graph is a graph whose edges are connected ~ Bern(p)
# This simulates a transition matrix for a random walk on an ER-p graph, where p = p_sparse.
RM_erdos <- function(N, p, stoch = T){
  P <- matrix(rep(NA, N * N), ncol = N)  # create [N x N] transition matrix
  for(i in 1:N){
    row <- runif(N,0,1) # Generate a uniform row of probabilites
    num_zeros <- rbinom(1,N,p) # Sample number of zeros so degree i ~ Bin(n,p)
    # Choose vertices to nullify based on sampled degree of vertex
    choices <- sample(1:N, num_zeros)
    row[choices] <- 0
    P[i,] <- row
    # Normalize rows
    for(i in 1:nrow(P)){
      row <- P[i, ]
      P[i,] <- row/sum(row)
    }
    # If the matrix is truly stochastic, rows with all zeros will have their diagonal become 1
    if(stoch){  
      # Set diagonal such that rows sum to 1
      diag <- rep(0, ncol(P))
      for(i in 1:nrow(P)){
        row <- P[i, ]
        diag[i] <- (1 - sum(nondiagonal_entries(row, i)))
      }
      diag(P) <- diag
    }
  }
  P # Return the matrix
}

#=================================================================================#
#                            STOCHASTIC ROW FUNCTIONS
#=================================================================================#

# generates stochastic rows of size N
r_stochastic <- function(N){
  prob <- runif(N,0,1)
  prob/sum(prob) # normalize
}

# generates same rows as in r_stochastic(N), but with introduced random sparsity
r_zeros <- function(N){
  prob <- runif(N,0,1)
  num_zeros <- sample(1:(N-1),1) # At most N-1 zeros, as to ensure stochastic property
  choices <- sample(1:N, num_zeros) # Choose edges to disconnect
  prob[choices] <- 0
  prob/sum(prob) # normalize
}

#=========================================================================#
#                             HELPER FUNCTIONS
#=========================================================================#

# Manually make equate the entries in the upper triangle to the conjugate of those in the lower triangle of the matrix
make_hermitian <- function(P){
  # Run over entry of the matrix
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      # Restrict view to one of the triangles (i < j): Lower Triangle
      if(i < j){P[i,j] <- Conj(P[j,i])} # Equalize lower and upper triangles, making conjugate if complex
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
