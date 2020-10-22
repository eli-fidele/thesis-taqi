#=================================================================================#
#                       RANDOM SYMMETRIC STOCHASTIC MATRICES
#=================================================================================#

rand_M_symm <- function(M,row_fn){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  # fill rows
  for(i in 1:M){
    P[i,] = row_fn(M,p)
  }
  P
}

#=================================================================================#
#                       RANDOM SPARSE STOCHASTIC MATRICES
#=================================================================================#

# initialize random P
rand_M_sparse <- function(M, p_sparse){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  for(i in 1:M){
    P[i,] = r_sparse(M,p)
  }
  P
}

# generates rows of size P which are valid probability distributions
r_sparse <- function(M,p){
  prob <- runif(M,0,1)
  num_zeros <- rbinom(1,M,p)
  choices <- sample(1:M, num_zeros)
  prob[choices] <- 0
  prob/sum(prob) # return normalized random row vector
}

#=================================================================================#
#                          RANDOM STOCHASTIC MATRICES
#=================================================================================#

# Possibly: instead of making seperate rand_M functions, feed list(args) to function, with default = empty set

# Generate random stochastic matrix of size M, with choice of row function {r_stochastic, r_zeros}
rand_M <- function(M,row_fn){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  # fill rows
  for(i in 1:M){
    P[i,] = row_fn(M,p)
    }
  P
}

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
