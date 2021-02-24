
#=================================================================================#
#                             RANDOM MATRIX ENSEMBLES
#=================================================================================#


#' @title Generate an ensemble of normal random matrices
#'
#' @description Given the same arguments as RM_norm, this function returns an ensemble of that particular class of matrix.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties (usually in spectral statistics) at the ensemble level.
#'
#' @inheritParams RM_norm
#' @param ... any default-valued parameters taken as arguments by the RM_ variant of this function
#' @param size the number of matrices to be simulated in the ensemble
#' @return An ensemble (list) of normal matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of standard normal 3x3 matrices of size 20
#' ensemble <- RME_norm(N = 3, size = 20)
#'
RME_norm <- function(N, ..., size){lapply(X = rep(N, size), FUN = RM_norm, ...)}


#' @title Generate an ensemble of random beta matrices
#'
#' @description Given the same arguments as RM_norm, this function returns an ensemble of that particular class of matrix.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties (usually in spectral statistics) at the ensemble level.
#'
#' @inheritParams RM_beta
#' @param ... any default-valued parameters taken as arguments by the RM_ variant of this function
#' @param size the number of matrices to be simulated in the ensemble
#' @return An ensemble (list) of beta matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of standard normal 3x3 matrices of size 20
#' # ensemble <-
#'
RME_beta <- function(N, ..., size){lapply(X = rep(N, size), FUN = RM_beta, ...)}


#' @title Generate an ensemble of stochastic matrices
#'
#' @description Given the same arguments as RM_stoch, this function returns an ensemble of that particular class of matrix.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties (usually in spectral statistics) at the ensemble level.
#'
#' @inheritParams RM_stoch
#' @param ... any default-valued parameters taken as arguments by the RM_ variant of this function
#' @param size the number of matrices to be simulated in the ensemble
#' @return An ensemble (list) of stochastic matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of standard normal 3x3 matrices of size 20
#' # ensemble <-
#'
RME_stoch <- function(N, ..., size){lapply(X = rep(N, size), FUN = RM_stoch, ...)}


#' @title Generate an ensemble of Erdos-Renyi transition matrices
#'
#' @description Given the same arguments as RM_norm, this function returns an ensemble of that particular class of matrix.
#'   While random matrices usually do not exude unique properties on their own, they do indeed have
#'   deterministic properties (usually in spectral statistics) at the ensemble level.
#'
#' @inheritParams RM_erdos
#' @param ... any default-valued parameters taken as arguments by the RM_ variant of this function
#' @param size the number of matrices to be simulated in the ensemble
#' @return An ensemble (list) of Erdos-Renyi transition matrices as specified by the matrix arguments.
#' @examples
#'
#' # Generate an ensemble of standard normal 3x3 matrices of size 20
#' # ensemble <-
#'
RME_erdos <- function(N, ..., size){lapply(X = rep(N, size), FUN = RM_erdos, ...)}

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
RM_norm <- function(N, mean = 0, sd = 1, symm = F, cplx = F, herm = F){
  # Create [N x N] matrix with normally distributed entries
  P <- matrix(rnorm(N^2, mean, sd), nrow = N)
  # Make symmetric if prompted
  if(symm || herm){P <- .make_hermitian(P)}
  # Returns a matrix with complex (and hermitian) entries if prompted
  if(cplx){
    if(herm){
      P <- P + .make_hermitian(1i * RM_norm(N, mean, sd))
    } else{
      P <- P + 1i * RM_norm(N, mean, sd, symm = F) # If not hermitian, recursively add a "real" instance of the imaginary components.
    }
  }
  P # Return the matrix
}

#' @title Generate a Gaussian/Hermite \eqn{\beta}-matrix
#'
#' @description Gaussian-\eqn{\beta} ensemble matrices are matrices with normal entries and beta real number components.
#'   Using Dumitriu's tridiagonal model, this function is an implementation of the generalized, but not necessarily invariant,
#'   beta ensembles for \eqn{\beta} > 0.
#'
#' @param N number of dimensions of the square matrix
#' @param beta the value of the beta parameter for the beta ensemble
#' @param cplx indicates whether the matrix should have complex entries.
#' @return A random Hermite beta matrix with any integer parameter beta
#' @examples
#' P <- RM_beta(N = 3, beta = 4)
#' P <- RM_beta(N = 10, beta = 17)
#'
RM_beta <- function(N, beta, cplx = F){
  # Set the diagonal as a N(0,2) distributed row.
  P <- diag(rnorm(N, mean = 0, sd = sqrt(2)))
  # Set the off-1 diagonals as chi squared variables with df(beta), as given in Dumitriu's model
  df_seq <- beta*(N - seq(1, N-1)) # Get degrees of freedom sequence for offdigonal
  P[row(P) - col(P) == 1] <- P[row(P) - col(P) == -1] <- sqrt(rchisq(N-1, df_seq)) # Generate tridiagonal
  # Add complex entries, if prompted
  #if(cplx){P <- P + .make_hermitian((1i * RM_beta(N, beta)))}
  P <- P/sqrt(2) # Rescale the entries by 1/sqrt(2)
  P # Return the matrix
}

#' Generate a tridiagonal matrix with normal entries
#'
#' @param N number of dimensions of the square matrix
#' @param symm indicates whether the matrix should be symmetric; equal to its transpose.
#' @return A random tridiagonal matrix with N(0,2) diagonal and N(0,1) band.
#' @examples
#' P <- RM_trid(N = 3)
#'
#' # Symmetric tridiagonal matrix
#' P <- RM_trid(N = 9, symm = TRUE)
#'
RM_trid <- function(N, symm = F){
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
    P <- .make_hermitian(P) # Make lower and upper triangles equal to each other's conjugate transpose
    diag(P) <- rep(0, N) # Nullify diagonal
    for(i in 1:N){P[i, ] <- P[i, ]/sum(P[i, ])} # Normalize rows
    # Set diagonal to the diff. between 1 and the non-diagonal entry sums such that rows sum to 1
    diag <- vector("numeric", N)
    for(i in 1:N){diag[i] <- (1 - sum(.nondiagonal_entries(row = P[i, ], row_index = i)))}
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
  degree_vertex <- rbinom(1, N, 1-p) # Sample number of zeros so that degree of row/vertex i ~ Bin(n,p)
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
  P # Return Hermitian Matrix
}

# Return the non-diagonal entries of row i
.nondiagonal_entries <- function(row, row_index){row[which(1:length(row) != row_index)]}

