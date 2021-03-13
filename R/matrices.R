
#=================================================================================#
#                             RANDOM MATRIX ENSEMBLES
#=================================================================================#

#' @title Generate an ensemble of normal random matrices
#'
#' @description Given the same arguments as RM_norm, this function returns an ensemble of random normal matrices.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties at the ensemble level in terms of their spectral statistics.
#'
#' @inheritParams RM_unif
#' @param ... any default-valued parameters taken as arguments by RM_norm()
#' @param size the size of the ensemble (i.e. number of matrices)
#' @return An ensemble (list) of normal matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of standard normal 3x3 matrices of size 20
#' ensemble <- RME_norm(N = 3, size = 20)
#'
RME_unif <- function(N, min, max, ..., size){lapply(X = rep(N, size), FUN = RM_unif, min, max, ...)}

#' @title Generate an ensemble of normal random matrices
#'
#' @description Given the same arguments as RM_norm, this function returns an ensemble of random normal matrices.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties at the ensemble level in terms of their spectral statistics.
#'
#' @inheritParams RM_norm
#' @param ... any default-valued parameters taken as arguments by RM_norm()
#' @param size the size of the ensemble (i.e. number of matrices)
#' @return An ensemble (list) of normal matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of standard normal 3x3 matrices of size 20
#' ensemble <- RME_norm(N = 3, size = 20)
#'
RME_norm <- function(N, mean = 0, sd = 1, ..., size){lapply(X = rep(N, size), FUN = RM_norm, mean, sd, ...)}


#' @title Generate an ensemble of random beta matrices
#'
#' @description Given the same arguments as RM_norm, this function returns an ensemble of that particular class of matrix.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties at the ensemble level in terms of their spectral statistics.
#'
#' @inheritParams RM_beta
#' @param size the size of the ensemble (i.e. number of matrices)
#' @return An ensemble (list) of beta matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of 10x10 beta matrices with beta = 4 of size 100.
#' ensemble <- RME_beta(N = 10, beta = 4, size = 100)
#'
RME_beta <- function(N, beta, size){lapply(X = rep(N, size), FUN = RM_beta, beta)}


#' @title Generate an ensemble of stochastic matrices
#'
#' @description Given the same arguments as RM_stoch, this function returns an ensemble of random stochastic matrices.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties at the ensemble level in terms of their spectral statistics.
#'
#' @inheritParams RM_stoch
#' @param ... pass any default-valued parameters taken as arguments by RM_stoch()
#' @param size the size of the ensemble (i.e. number of matrices)
#' @return An ensemble (list) of stochastic matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of random 5x5 transition matrices of size 20.
#' ensemble <- RME_stoch(N = 5, size = 20)
#'
#' # Generate an ensemble of symmetric random 5x5 transition matrices of size 20.
#' ensemble <- RME_stoch(N = 5, symm = TRUE, size = 20)
#'
RME_stoch <- function(N, ..., size){lapply(X = rep(N, size), FUN = RM_stoch, ...)}


#' @title Generate an ensemble of Erdos-Renyi transition matrices
#'
#' @description Given the same arguments as RM_norm, this function returns an ensemble of random Erdos-Renyi stochastic matrices.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties at the ensemble level in terms of their spectral statistics.
#'
#' @inheritParams RM_erdos
#' @param ... any default-valued parameters taken as arguments by RM_erdos()
#' @param size the size of the ensemble (i.e. number of matrices)
#' @return An ensemble (list) of Erdos-Renyi transition matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of 10x10 Erdos-Renyi transition matrices of size 50 with p = 0.7
#' # ensemble <- RME_erdos(N = 10, p = 0.7, size = 50)
#'
RME_erdos <- function(N, p, ..., size){lapply(X = rep(N, size), FUN = RM_erdos, p, ...)}

#=================================================================================#
#                                RANDOM MATRICES
#=================================================================================#

#' @title Generate a uniform random matrix
#'
#' @description Uniform random matrices are matrices with uniformly distributed entries. They are an elementary type of random matrix.
#'
#' @param N number of dimensions of the square matrix
#' @param min minimum of the uniform distribution to be sampled from
#' @param max maximum of the uniform distribution to be sampled from
#' @param symm indicates whether the matrix should be symmetric (equal to its transpose).
#' @param cplx indicates whether the matrix should have complex entries.
#' @param herm indicates whether the matrix should be hermitian (equal to its conjugate transpose).
#'   Reserved for when complex = T, otherwise use symm = T.
#' @return A random matrix with uniformly distributed entries.
#' @examples
#' # Unif(1,2) distributed matrix
#' P <- RM_unif(N = 3, min = 1, max = 2)
#'
#' # Unif(0,5) distributed matrix with real symmetric entries
#' P <- RM_unif(N = 7, min = 0, max = 5, symm = TRUE)
#'
#' # Unif(0,1) distributed matrix with complex entries
#' Q <- RM_unif(N = 7, min = 0, max = 1, cplx = TRUE)
#'
#' # Unif(2,10) distributed matrix with hermitian complex entries
#' Q <- RM_unif(N = 5, min = 2, max = 10, cplx = TRUE, herm = TRUE)
#'
RM_unif <- function(N, min, max, symm = FALSE, cplx = FALSE, herm = FALSE){
  # Create [N x N] matrix with uniformly distributed entries
  P <- matrix(runif(N^2, min, max), nrow = N)
  # Make symmetric/hermitian if prompted
  if(symm || herm){P <- .makeHermitian(P)}
  # Returns a matrix with complex (and hermitian) entries if prompted
  if(cplx){
    Im_P <- (1i * RM_unif(N, min, max)) # Recursively add imaginary components as 1i * instance of real-valued matrix.
    if(herm){P <- P + .makeHermitian(Im_P)}  # Make imaginary part hermitian if prompted
    else{P <- P + Im_P}
  }
  P # Return the matrix
}

#=================================================================================#
#                           NORMAL RANDOM MATRICES
#=================================================================================#

#' @title Generate a normal random matrix
#'
#' @description Normal random matrices are matrices with normally distributed entries. These matrices
#'  are extensively studied in random matrix theory.
#'
#' @param N number of dimensions of the square matrix
#' @param mean mean of the normal distribution of entries
#' @param sd standard deviation of the normal distribution of entries
#' @param symm indicates whether the matrix should be symmetric (equal to its transpose).
#'   Reserved for when complex = F, otherwise use hermitian = T.
#' @param cplx indicates whether the matrix should have complex entries.
#' @param herm indicates whether the matrix should be hermitian (equal to its conjugate transpose).
#'   Reserved for when complex = T, otherwise use symm = T.
#' @return A random matrix with normally distributed entries.
#' @examples
#' # N(1,2) distributed matrix
#' P <- RM_norm(N = 3, mean = 1, sd = 2)
#'
#' # N(0,5) distributed matrix with real symmetric entries
#' P <- RM_norm(N = 7, sd = 5, symm = TRUE)
#'
#' # 7x7 standard normal matrix with complex entries
#' Q <- RM_norm(N = 7, cplx = TRUE)
#'
#' # N(2,1) distributed matrix with hermitian complex entries
#' Q <- RM_norm(N = 5, mean = 2, cplx = TRUE, herm = TRUE)
#'
RM_norm <- function(N, mean = 0, sd = 1, symm = FALSE, cplx = FALSE, herm = FALSE){
  # Create [N x N] matrix with normally distributed entries
  P <- matrix(rnorm(N^2, mean, sd), nrow = N)
  # Make symmetric/hermitian if prompted
  if(symm || herm){P <- .makeHermitian(P)}
  # Returns a matrix with complex (and hermitian) entries if prompted
  if(cplx){
    Im_P <- (1i * RM_norm(N, mean, sd)) # Recursively add imaginary components as 1i * instance of real-valued matrix.
    if(herm){P <- P + .makeHermitian(Im_P)}  # Make imaginary part hermitian if prompted
    else{P <- P + Im_P}
  }
  P # Return the matrix
}

#' @title Generate a Hermite \eqn{\beta}-matrix
#'
#' @description Hermite-\eqn{\beta} ensemble matrices are matrices with normal entries and beta real number components.
#'   Using Dumitriu's tridiagonal model, this function is an implementation of the generalized, but not necessarily invariant,
#'   beta ensembles for \eqn{\beta} > 0.
#'
#' @param N number of dimensions of the square matrix
#' @param beta the value of the beta parameter for the beta ensemble
#' @return A random Hermite beta matrix with any integer parameter beta
#' @examples
#'
#' # Generate a 3x3 random beta matrix with beta = 4
#' P <- RM_beta(N = 3, beta = 4)
#'
#' # Generate a 10x10 random beta matrix with beta = 21
#' P <- RM_beta(N = 10, beta = 21)
#'
RM_beta <- function(N, beta){
  # Set the diagonal as a N(0,2) distributed row.
  P <- diag(rnorm(N, mean = 0, sd = sqrt(2)))
  # Set the off-1 diagonals as chi squared variables with df(beta), as given in Dumitriu's model
  df_seq <- beta*(N - seq(1, N-1)) # Get degrees of freedom sequence for offdigonal
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- sqrt(rchisq(N-1, df_seq)) # Generate tridiagonal
  P <- P/sqrt(2) # Rescale the entries by 1/sqrt(2)
  P # Return the matrix
}

#' Generate a tridiagonal matrix with normal entries
#'
#' @param N number of dimensions of the square matrix
#' @param symm indicates whether the matrix should be symmetric; equal to its transpose.
#' @return A random tridiagonal matrix with N(0,2) diagonal and N(0,1) band.
#' @examples
#'
#' # Generate a 3x3 standard normal tridiagonal matrix
#' P <- RM_trid(N = 3)
#'
#' # Symmetric tridiagonal matrix
#' P <- RM_trid(N = 9, symm = TRUE)
#'
RM_trid <- function(N, symm = FALSE){
  diagonal <- rnorm(n = N, 0, 2)
  P <- diag(diagonal)
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- rnorm(n = N, 0, 1)
  P # Return the matrix
}

#=================================================================================#
#                           STOCHASTIC RANDOM MATRICES
#=================================================================================#

#' @title Generate a random stochastic matrix
#'
#' @description A (row-)stochastic matrix is a matrix whose rums sum to 1. There is a natural one-to-one corrospondence between
#'   stochastic matrices and Markov Chains; this is so when its i,j entry represent the transition probability from state i to state j.
#'
#' @param N number of dimensions of the square matrix
#' @param symm indicates whether the matrix should be symmetric; equal to its transpose.
#' @param sparsity indicates whether the matrix should add some arbitrary sparsity (zeros)
#' @return A random stochastic matrix.
#' @examples
#' P <- RM_stoch(N = 3)
#' P <- RM_stoch(N = 9, sparsity = TRUE)
#' Q <- RM_stoch(N = 9, symm = TRUE)
#' Q <- RM_stoch(N = 9, symm = TRUE, sparsity = TRUE)
#'
RM_stoch <- function(N, symm = F, sparsity = F){
  if(sparsity){row_fxn <- .stoch_row_zeros} else {row_fxn <- .stoch_row} # Choose row function
  # Generate the [N x N] stochastic matrix stacking N stochastic rows (using the chosen function)
  P <- do.call("rbind", lapply(X = rep(N, N), FUN = row_fxn))
  if(symm){ # Make symmetric (if prompted)
    P <- .makeHermitian(P) # Make lower and upper triangles equal to each other's conjugate transpose
    diag(P) <- rep(0, N) # Nullify diagonal
    for(i in 1:N){P[i, ] <- P[i, ]/sum(P[i, ])} # Normalize rows
    # Set diagonal to the diff. between 1 and the non-diagonal entry sums such that rows sum to 1
    diag <- vector("numeric", N)
    for(i in 1:N){diag[i] <- (1 - sum(.offdiagonalEntries(row = P[i, ], row_index = i)))}
    diag(P) <- diag
  }
  P # Return the matrix
}

#' @title Generate a random stochastic matrix for a walk on an Erdos-Renyi graph
#'
#' @description An Erdos-Renyi Graph is a graph whose edges are connected ~ Bern(p).
#'  Hence, its transition matrix will have nonzero entries with that probability.
#'  So, we can alternatively think of the transition matrix for such walk as a stochastic matrix with parameterized sparsity.
#'
#' @param N number of dimensions of the square matrix
#' @param p the probability two vertices are connected in an Erdos-Renyi graph.
#' @param stoch indicates whether the matrix should be stochastic;
#'   changing this parameter may lead to the function returning invalid transition matrices.
#' @return A random stochastic matrix corrosponding to a walk on an Erdos-Renyi graph with probability p.
#' @examples
#' # Very sparse graph
#' P <- RM_erdos(N = 3, p = 0.2)
#'
#' P <- RM_erdos(N = 9, p = 0.6)
#'
#' # Completely connected graph
#' P <- RM_erdos(N = 5, p = 1)
#'
#' # May yield invalid transition matrix.
#' Q <- RM_erdos(N = 5, p = 0.3, stoch = FALSE)
#'
RM_erdos <- function(N, p, stoch = T){
  # Generate an [N x N] Erdos-Renyi walk stochastic matrix by stacking N p-stochastic rows
  P <- do.call("rbind", lapply(X = rep(N, N), FUN = .stoch_row_erdos, p = p))
  # If the matrix is to be truly stochastic, map rows with all zeros to have diagonal entry 1
  if(stoch){
    # Set diagonal to ensure that rows sum to 1
    diag <- rep(0, N)
    for(i in 1:N){diag[i] <- (1 - sum(.offdiagonalEntries(row = P[i, ], row_index = i)))}
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
  degree_vertex <- rbinom(1, N, 1-p) # Sample number of zeros so that degree of row/vertex i ~ Bin(n,p)
  row[sample(1:N, degree_vertex)] <- 0 # Choose edges to sever and sever them
  if(sum(row) != 0){row/sum(row)} else{row} # Return normalized row only if non-zero (cannot divide by 0)
}


#=========================================================================#
#                             HELPER FUNCTIONS
#=========================================================================#

# Manually make equate the entries in the upper triangle to the conjugate of those in the lower triangle of the matrix
.makeHermitian <- function(P){
  # Run over entry of the matrix
  for(i in 1:nrow(P)){
    for(j in 1:ncol(P)){
      # Restrict view to one of the triangles (i < j): Lower Triangle
      if(i < j){P[i,j] <- Conj(P[j,i])} # Equalize lower and upper triangles, making conjugate if complex
    }
  }
  P # Return Hermitian Matrix
}

# Return the off-diagonal entries of row i 
.offdiagonalEntries <- function(row, row_index){row[which(1:length(row) != row_index)]}

