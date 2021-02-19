
#=================================================================================#
#                           NORMAL RANDOM MATRICES 
#=================================================================================#

RM_norm <- function(N, mean = 0, sd = 1, symm = F, complex = F, hermitian = F){
  # Create [n x n] matrix with normally distributed entries
  P <- matrix(rnorm(N^2, mean, sd), nrow = N)  
  # Make symmetric if prompted
  if(symm || hermitian){P <- .make_hermitian(P)}
  # Returns a matrix with complex entries if prompted
  if(complex){
    if(hermitian){
      P <- P + .make_hermitian(1i * RM_normal(N, mean, sd, symm = T))
    } else{
      P <- P + 1i * RM_normal(N, mean, sd, symm = F)
    }
  }
  P <- P/sqrt(N) # Rescale the matrix
  P # Return the matrix
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
  if(complex){P <- P + .make_hermitian((1i * RM_beta(N, beta)))}
  P <- P/sqrt(2) # Rescale the entries by 1/sqrt(2)
  P # Return the matrix
}

# Generate a tridiagonal matrix with normal entries
RM_trid <- function(N, symm = F){
  diagonal <- rnorm(N = N, 0, 2)
  P <- diag(diagonal)
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- rnorm(n = n, 0, 1)
  P # Return the matrix
}

#=================================================================================#
#                           STOCHASTIC RANDOM MATRICES 
#=================================================================================#

# Generate random stochastic matrix of size n, with choice of row function {r_stoch, r_stoch_zeros}
RM_stoch <- function(N, symm = F, sparsity = F){
  if(sparsity){row_fxn <- .stoch_row_zeros} else {row_fxn <- .stoch_row} # Choose row function
  # Generate the [N x N] stochastic matrix stacking N stochastic rows (using the chosen function)
  P <- do.call("rbind", lapply(X = rep(N, N), FUN = row_fxn))
  if(symm){ # Make symmetric (if prompted)
    P <- .make_hermitian(P) # Make lower and upper triangles equal
    diag(P) <- rep(0, N) # Nullify diagonal
    for(i in 1:N){P[i, ] <- P[i, ]/sum(P[i, ])} # Normalize rows
    # Set diagonal to the diff. between 1 and the non-diagonal entry sums such that rows sum to 1
    diag <- vector("numeric", N)
    for(i in 1:N){diag[i] <- (1 - sum(.nondiagonal_entries(row = P[i, ], row_index = i)))}
    diag(P) <- diag
  }
  P # Return the matrix
}

# An Erdos-Renyi Graph is a graph whose edges are connected ~ Bern(p).
# This simulates a transition matrix for a random walk on an ER-p graph, where p = p_sparse.
RM_erdos <- function(N, p, stoch = T){
  # Generate an [N x N] Erdos-Renyi walk stochastic matrix by stacking N p-stochastic rows (using the chosen function)
  P <- do.call("rbind", lapply(X = rep(N, N), FUN = .stoch_row_erdos, p = p))
  # If the matrix is to be truly stochastic, map rows with all zeros to have diagonal entry 1
  if(stoch){  
    # Set diagonal to ensure that rows sum to 1
    diag <- rep(0, N)
    for(i in 1:N){diag[i] <- (1 - sum(.nondiagonal_entries(row = P[i, ], row_index = i)))}
    diag(P) <- diag
  }
  P # Return the matrix
}

#=================================================================================#
#                            STOCHASTIC ROW FUNCTIONS
#=================================================================================#

# Generates stochastic rows of size N
.stoch_row <- function(N){
  row <- runif(N,0,1) # Sample probability distribution
  row/sum(row) # Return normalized row
}

# Generates same rows as in r_stoch(N), but with introduced random sparsity
.stoch_row_zeros <- function(N){
  row <- runif(N,0,1)
  degree_vertex <- sample(1:(N-1), size = 1) # Sample a degree of at least 1, as to ensure row is stochastic
  row[sample(1:N, size = degree_vertex)] <- 0 # Choose edges to sever and sever them
  row/sum(row) # Return normalized row
}

# Generates a stochastic row with parameterized sparsity of p
.stoch_row_erdos <- function(N, p){
  row <- runif(N,0,1) # Generate a uniform row of probabilites
  degree_vertex <- rbinom(1,N,1-p) # Sample number of zeros so that degree of row/vertex i ~ Bin(n,p)
  row[sample(1:N, degree_vertex)] <- 0 # Choose edges to sever and sever them
  if(sum(row) != 0){row/sum(row)} else{row} # Return normalized row only if non-zero (cannot divide by 0)
}


#=========================================================================#
#                             HELPER FUNCTIONS
#=========================================================================#

# Manually make equate the entries in the upper triangle to the conjugate of those in the lower triangle of the matrix
.make_hermitian <- function(P){
  # Run over entry of the matrix
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      # Restrict view to one of the triangles (i < j): Lower Triangle
      if(i < j){P[i,j] <- Conj(P[j,i])} # Equalize lower and upper triangles, making conjugate if complex
    }
  }
  P # Return Symmetric Matrix
}

# Return the non-diagonal entries of row i
.nondiagonal_entries <- function(row, row_index){
  row[which(1:length(row) != row_index)]
  }

# Check if a matrix is stochastic
.isStochastic <- function(P){
  row_is_stoch <- rep(F, nrow(P))
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    row_is_stoch[i] <- (row_sum == 1)
  }
  !(F %in% row_is_stoch)
}

# Returns proportion of positive entries of any matrix P
.pos_entries <- function(P, zero = F){
  if(!zero){pos_entries <- length(matrix(P[P[,] > 0], nrow = 1))}
  else{pos_entries <- length(matrix(P[P[,] >= 0], nrow = 1))}
  pos_entries/(length(P))   
}
