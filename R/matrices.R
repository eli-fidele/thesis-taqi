
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
  if(symm == T){
    P[lower.tri(P)] <- P[upper.tri(P)]
    }
  # Return the matrix
  P
}

# Generate random stochastic matrix of size M, with choice of row function {r_stochastic, r_zeros}
RM_stoch <- function(M, symm = F, sparsity = F){
  P <- matrix(rep(NA, M * M), ncol = M)  # create [M x M] transition matrix
  if(sparsity){row_fn <- r_zeros} else {row_fn <- r_stochastic} # choose row function
  # Generate rows
  for(i in 1:M){P[i,] <- row_fn(M)}
  # Make symmetric (if prompted)
  if(symm == T){
    # Add transpose to make it symmetric
    P[lower.tri(P)] <- P[upper.tri(P)]
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
#                         RANDOM MATRIX DIAGNOSTICS 
#=================================================================================#

nondiagonal_entries <- function(row, row_index){
  indices <- data.frame(idx = 1:length(row))
  indices <- indices %>% filter(idx != row_index)
  # return the row with the given indices
  row[as.numeric(indices[,])]
}

# returns proportion of positive entries of any matrix P
pos_entries <- function(P){
  pos_entries <- length(matrix(P[P[,] > 0], nrow = 1))
  pos_entries/(length(P))   
}

# Vectorize matrix entries to study their distribution
vectorize_matrix <- function(P){as.vector(P)}

# Check if a matrix is symmetric 
is_symmetric <- function(P){isSymmetric(P)}

# Check if a matrix is stochastic
is_row_stochastic <- function(P){
  row_is_stoch <- rep(F, nrow(P))
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    row_is_stoch[i] <- (row_sum == 1)
  }
  !(F %in% row_is_stoch)
}


