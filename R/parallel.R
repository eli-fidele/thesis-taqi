
#=================================================================================#
#                           EIGENVALUE DISPERSION (PARALLEL)
#=================================================================================#
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
#                           EIGENVALUE DISPERSION (PARALLEL)
#=================================================================================#
#' @title Obtain the eigenvalue spacings of a matrix or ensemble of matrices.
#'
#' @description Returns a vector of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#' @examples
#' 
#' # Eigenvalue dispersion in parallel
#' P <- RM_norm(N = 100, size = 500)
#' #disp_P <- dispersion_parallel(P)
dispersion_parallel <- function(array, pairs = NA, sort_norms = TRUE, singular = FALSE, norm_pow = 1){ #sortNorms? orderByNorms? pair_scheme?
  digits <- 4 # Digits to round values to
  pairs <- .parsePairs(pairs, array) # Parse input and generate pair scheme (default NA), passing on array for dimension and array type inference
  # Array is a matrix; call function returning dispersion for singleton matrix
  if(class(array) == "matrix"){
    .dispersion_matrix(array, pairs, sort_norms, singular, norm_pow, digits)
  }
  # Array is an ensemble; recursively row binding each matrix's dispersions
  else if(class(array) == "list"){
    furrr::future_map_dfr(array, .dispersion_matrix, pairs, sort_norms, singular, norm_pow, digits)
  }
}