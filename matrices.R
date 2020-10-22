#=================================================================================#
#                         SYMMETRIC TRIDIAGONAL MATRICES
#=================================================================================#

rand_M_trid <- function(M){
  diagonal <- rnorm(n = M, 0, 2)
  P <- diag(diagonal)
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- rnorm(n = M, 0, 1)
  P
}

#=================================================================================#
#                           NORMAL SYMMETRIC MATRICES
#=================================================================================#

rand_M_symm <- function(M, mu, sd){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  # fill rows
  for(i in 1:M){
    P[i,] = r_normal(M, mu, sd)
  }
  P %*% t(P) 
}

r_normal <- function(M, mu, sd){
  rnorm(n = M, mean = mu, sd = sd)
}

#=================================================================================#
#                          SYMMETRIC STOCHASTIC MATRICES
#=================================================================================#

rand_M_symm_stoch <- function(M,row_fn){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  # fill rows
  for(i in 1:M){
    P[i,] = row_fn(M)
  }
  P %*% t(P) 
}

#=================================================================================#
#               p-SPARSE STOCHASTIC MATRICES (ERDOS-RENYI GRAPHS)
#=================================================================================#

# initialize random P
rand_M_erdos <- function(M, p_sparse){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  for(i in 1:M){
    P[i,] = r_sparse(M,p_sparse) # p_sparse is a probability between [0,1) 
  }                              # so edges are connected ~ Bern(p) 
  P
}

# generates rows of size P which are valid probability distributions
r_sparse <- function(M,p){
  prob <- runif(M,0,1)
  num_zeros <- rbinom(1,M,p) # Isomorphic to Erdos-Renyi graphs!
  choices <- sample(1:M, num_zeros)
  prob[choices] <- 0
  prob/sum(prob) # return normalized random row vector
}

#=================================================================================#
#                             STOCHASTIC MATRICES
#=================================================================================#

# Possibly: instead of making seperate rand_M functions, feed list(args) to function, with default = empty set

# Generate random stochastic matrix of size M, with choice of row function {r_stochastic, r_zeros}
rand_M_stoch <- function(M, row_fn){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  # fill rows
  for(i in 1:M){
    P[i,] = row_fn(M)
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
