
#=================================================================================#
#                             RANDOM MATRIX FUNCTIONS
#=================================================================================#

RM_normal <- function(M, normal_args = c(0,1), symm = F){
  mu <- normal_args[1]
  sd <- normal_args[2]
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  # fill rows
  for(i in 1:M){
    P[i,] <- rnorm(n = M, mean = mu, sd = sd)
  }
  if(symm == T){P <- P %*% t(P)}
  # Return the matrix
  P
}

# Generate random stochastic matrix of size M, with choice of row function {r_stochastic, r_zeros}
RM_stoch <- function(M, symm = F, sparsity = F){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  if(sparsity){row_fn <- r_zeros} else {row_fn <- r_stochastic}
  # fill rows
  for(i in 1:M){
    P[i,] <- row_fn(M)
  }
  if(symm == T){P <- P %*% t(P)}
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
    if (stoch == T){curr_row <- curr_row/sum(curr_row)} # Normalize if to be stochastic
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
#                         RANDOM MATRIX HELPER FUNCTIONS
#=================================================================================#

# Vectorize matrix entries to study their distribution

vectorize_matrix <- function(P){as.vector(P)}



