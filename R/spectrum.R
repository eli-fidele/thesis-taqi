
# This script includes functions that help extract and plot eigenvalues of matrices.

#=================================================================================#
#                           EIGENVALUE SPECTRUM (PARALLEL)
#=================================================================================#
#' 
#' @title Obtain the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble.
#'
#' @inheritParams spectrum
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#' @examples
#' 
#' # Eigenvalue spectrum computed in parallel
#' P <- RME_norm(N = 100, size = 500)
#' #spectrum_P <- spectrum_parallel(P)
#' 
spectrum_parallel <- function(array, components = TRUE, sort_norms = TRUE, singular = FALSE, order = NA){
  digits <- 4 # Digits to round values to
  # Array is a matrix; call function returning eigenvalues for singleton matrix
  if(class(array) == "matrix"){
    .spectrum_matrix(array, components, sort_norms, singular, order, digits)
  }
  # Array is an ensemble; recursively row binding each matrix's eigenvalues
  else if(class(array) == "list"){
    furrr::future_map_dfr(array, .spectrum_matrix, components, sort_norms, singular, order, digits)
  }
}

#=================================================================================#
#                              EIGENVALUE SPECTRUM 
#=================================================================================#

#' @title Obtain the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be returned
#' @param components returns the array with resolved real and imaginary components if TRUE; otherwise returns complex-valued eigenvalues
#' @param sort_norms sorts the eigenvalue spectrum by its norms, otherwise sorts by sign
#' @param singular get the singular values of the matrix (i.e. square root of the eigenvalues of the matrix times its transpose)
#' @param order get eigenvalues with that given order (norm ranking); order 1 represents largest, order N represents smallest (where N is the number of eigenvalues).
#'   If uninitialized, returns the entire spectrum.
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#' @examples
#'
#' # Eigenvalue spectrum of a random normal matrix
#' P <- RM_norm(N = 5)
#' spectrum_P <- spectrum(P)
#'
#' Q <- matrix(runif(2^2), ncol = 2)
#' spectrum_Q <- spectrum(Q)
#'
#' # Eigenvalue spectra of ensemble matrices
#' ensemble <- RME_norm(N = 3, size = 10)
#' ensemble_spectrum <- spectrum(ensemble)
#'
spectrum <- function(array, components = TRUE, sort_norms = TRUE, singular = FALSE, order = NA){
  digits <- 4 # Digits to round values to
  # Array is an ensemble; recursively row binding each matrix's eigenvalues
  if(class(array) == "list"){
    purrr::map_dfr(array, .spectrum_matrix, components, sort_norms, singular, order, digits)
  }
  # Array is a matrix; call function returning eigenvalues for singleton matrix
  else{
    .spectrum_matrix(array, components, sort_norms, singular, order, digits)
  }
}

# Helper function returning tidied eigenvalue array for a matrix
.spectrum_matrix <- function(P, components, sort_norms, singular, order, digits = 4){
  # If prompted for singular values, then take the product of the matrix and its tranpose instead
  if(singular){P <- P %*% t(P)}
  # Get the sorted eigenvalue spectrum of the matrix
  eigenvalues <- eigen(P)$values # Compute the eigenvalues of P
  if(singular){eigenvalues <- sqrt(eigenvalues)} # Take the square root of the eigenvalues
  if(sort_norms){eigenvalues <- .sort_by_norm(eigenvalues)} # Order the eigenvalue spectrum by norm rather than sign
  else{eigenvalues <- sort(eigenvalues, decreasing = TRUE)} # Else, sort by sign.
  # If uninitialized, get eigenvalues of all orders; Otherwise, concatenate so single inputs become vectors
  if(class(order) == "logical"){order <- 1:nrow(P)} else{order <- c(order)}
  purrr::map_dfr(order, .resolve_eigenvalue, eigenvalues, components, digits) # Get the eigenvalues
}

# Read and parse an eigenvalue from an eigen(P)$value array
.resolve_eigenvalue <- function(order, eigenvalues, components, digits){
  eigenvalue <- eigenvalues[order] # Read from eigen(P)$values
  # Get norm and order columns (will unconditionally be returned)
  norm_and_order <- data.frame(Norm = abs(eigenvalue), Order = order)
  # If components are requested, resolve parts into seperate columns and cbind to norm and order
  if(components){evalue <- cbind(data.frame(Re = Re(eigenvalue), Im = Im(eigenvalue)), norm_and_order)}
  else{evalue <- cbind(data.frame(Eigenvalue = eigenvalue), norm_and_order)}
  evalue <- round(evalue, digits) # Round entries
  evalue # Return resolved eigenvalue
}

#=================================================================================#
#                              ORDER SORTING SCHEMES
#=================================================================================#

# Sort an array of numbers by their norm (written for eigenvalue sorting)
.sort_by_norm <- function(eigenvalues){
  (data.frame(eigenvalue = eigenvalues, norm = abs(eigenvalues)) %>% arrange(desc(norm)))$eigenvalue
  }


