
#=================================================================================#
#                           EIGENVALUE SPECTRUM (PARALLEL)
#=================================================================================#

#' @title Obtain the ordered eigenvalue spectrum of a matrix or ensemble of matrices concurrently (faster).
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble concurrently (faster).
#'
#' @inheritParams spectrum
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#'
#' @examples
#' # Generate a random matrix
#' P <- RME_norm(N = 20, size = 100)
#' # Compute the spectrum concurrently (faster)
#' #spec_P <- spectrum_par(P)
#'
spectrum_par <- function(array, norm_order = TRUE, singular = FALSE, components = TRUE, order = NA){
  # Digits to round values to
  digits <- 4
  # Set up futures
  future::plan(future::multisession)
  # Get the type of array
  array_class <- .arrayClass(array)
  # Compute the spectrum
  if(array_class == "ensemble"){
    # Array is an ensemble; recursively row bind each matrix's eigenvalues
    furrr::future_map_dfr(array, .spectrum_matrix, norm_order, singular, components, order, digits)
  }
  else if(array_class == "matrix"){
    # Array is a matrix; call function returning eigenvalues for a singleton matrix
    .spectrum_matrix(array, norm_order, singular, components, order, digits)
  }
}

#=================================================================================#
#                               DISPERSION (PARALLEL)
#=================================================================================#

#' @title Obtain the eigenvalue spacings of a matrix or ensemble of matrices.
#' @description Returns a vector of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#'
#' @examples
#' # Eigenvalue dispersion in parallel
#' P <- RME_norm(N = 20, size = 100)
#' #disp_P <- dispersion_par(P)
#'
dispersion_par <- function(array, pairs = NA, norm_order = TRUE, singular = FALSE, pow_norm = 1){
  # Digits to round values to
  digits <- 4
  # Set up futures
  future::plan(future::multisession)
  # Get the type of array
  array_class <- .arrayClass(array)
  # Parse input and generate pair scheme (default NA), passing on array for dimension and array type inference
  pairs <- .parsePairs(pairs, array, array_class)
  # Compute the dispersion
  if(array_class == "ensemble"){
    # Array is an ensemble; recursively row binding each matrix's dispersions
    furrr::future_map_dfr(array, .dispersion_matrix, pairs, norm_order, singular, pow_norm, digits)
  }
  else if(array_class == "matrix"){
    # Array is a matrix; call function returning dispersion for singleton matrix
    .dispersion_matrix(array, pairs, norm_order, singular, pow_norm, digits)
  }
}