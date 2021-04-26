
# spectrum.R : A script file including functions to help compute the spectra
# of random matrices and random matrix ensembles.

#=================================================================================#
#                              EIGENVALUE SPECTRUM
#=================================================================================#

#' @title Obtain the ordered eigenvalue spectrum of a matrix or ensemble of matrices.
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be returned
#' @param norm_order sorts the eigenvalue spectrum by its norms if TRUE, otherwise sorts them by sign
#' @param singular return the singular values of the matrix or matrix ensemble
#' @param components returns the array with resolved real and imaginary components if TRUE, otherwise returns complex-valued eigenvalues
#' @param order an integer or integer vector of which eigenvalue orders to return; order 1 representing the largest, order N represents smallest (where N is the number of eigenvalues).
#'   If uninitialized, defaults to returning the entire spectrum.
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#'
#' @examples
#' # Eigenvalue spectrum of a random normal matrix
#' P <- RM_norm(N = 5)
#' spec_P <- spectrum(P)
#'
#' Q <- matrix(runif(2^2), ncol = 2)
#' spec_Q <- spectrum(Q)
#'
#' # Eigenvalue spectra of ensemble matrices
#' ens <- RME_norm(N = 3, size = 10)
#' spec_ens <- spectrum(ens)
#'
spectrum <- function(array, norm_order = TRUE, singular = FALSE, components = TRUE, order = NA){
  # Digits to round values to
  digits <- 4
  # Get the type of array
  array_class <- .arrayClass(array)
  # Depending on the array class, call the appopriate functions
  if(array_class == "ensemble"){
    # For ensembles, iteratively rbind() each matrix's spectrum
    purrr::map_dfr(array, .spectrum_matrix, norm_order, singular, components, order, digits)
  }
  else if(array_class == "matrix"){
    # For matrices, call the function returning the ordered spectrum for a singleton matrix
    .spectrum_matrix(array, norm_order, singular, components, order, digits)
  }
}

#=================================================================================#
# Helper function returning tidied eigenvalue array for a matrix
.spectrum_matrix <- function(P, norm_order, singular, components, order, digits = 4){
  # For singular values, take P as product of the itself and its tranpose
  if(singular){P <- P %*% t(P)}
  # Get the eigenvalues of P
  eigenvalues <- eigen(P, only.values = TRUE)$values
  # Take the square root of the eigenvalues to obtain singular values
  if(singular){eigenvalues <- sqrt(eigenvalues)}
  # Sort the eigenvalues to make it an ordered spectrum
  eigenvalues <- .sortValues(eigenvalues, norm_order)
  # If uninitialized, get eigenvalues of all orders; otherwise, use c() so singletons => vectors
  if(class(order) == "logical"){ order <- 1:nrow(P) } else{ order <- c(order) }
  # Return the spectrum of the matrix
  purrr::map_dfr(order, .resolve_eigenvalue, eigenvalues, components, digits)
}

#=================================================================================#
# Read and parse an eigenvalue from a sorted eigenvalue array
.resolve_eigenvalue <- function(order, eigenvalues, components, digits){
  # Read from a sorted eigenvalue array at that order
  eigenvalue <- eigenvalues[order]
  # Get norm and order columns
  features <- data.frame(Norm = abs(eigenvalue), Order = order)
  if(components){
    # If components are sought, resolve the eigenvalue into seperate columns first
    res <- cbind(data.frame(Re = Re(eigenvalue), Im = Im(eigenvalue)), features)
  } else{
    # Otherwise, don't resolve the eigenvalue components
    res <- cbind(data.frame(Eigenvalue = eigenvalue), features)
  }
  # Round entries and return the resolved eigenvalue
  res <- round(res, digits)
  return(res)
}

#=================================================================================#
#                           SPECTRUM: HELPER FUNCTIONS
#=================================================================================#

# Parses an array to see classify it as a matrix or an ensemble of matrices.
.arrayClass <- function(array){
  # Sample an element from the array and get its class
  elem <- array[[1]]
  types <- class(elem)
  # Classify it by analyzing the element class
  if("numeric" %in% types || "complex" %in% types){
    return("matrix")
  }
  else if("matrix" %in% types){
    return("ensemble")
  }
}

# Sort an array of numbers by their norm (written for eigenvalue sorting)
.sortValues <- function(vals, norm_order){
  values <- data.frame(value = vals)
  # If asked to sort by norms, arrange by norm and return
  if(norm_order){
    values$norm <- abs(values$value)
    values <- values %>% arrange(desc(norm))
    # Return the norm-sorted values
    values$value
  }
  # Otherwise, sort by sign and return
  else{ sort(vals, decreasing = TRUE) }
}
