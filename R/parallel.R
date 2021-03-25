
spectrum <- function(array, components = TRUE, sort_norms = TRUE, singular = FALSE, order = NA){
  digits <- 4 # Digits to round values to
  # Array is a matrix; call function returning eigenvalues for singleton matrix
  if(class(array) == "matrix"){
    .spectrum_matrix(array, components, sort_norms, singular, order, digits)
  }
  # Array is an ensemble; recursively row binding each matrix's eigenvalues
  else if(class(array) == "list"){
    purrr::map_dfr(array, .spectrum_matrix, components, sort_norms, singular, order, digits)
  }
}

dispersion <- function(array, pairs = NA, sort_norms = TRUE, singular = FALSE, norm_pow = 1){ #sortNorms? orderByNorms? pair_scheme?
  digits <- 4 # Digits to round values to
  pairs <- .parsePairs(pairs, array) # Parse input and generate pair scheme (default NA), passing on array for dimension and array type inference
  # Array is a matrix; call function returning dispersion for singleton matrix
  if(class(array) == "matrix"){
    .dispersion_matrix(array, pairs, sort_norms, singular, norm_pow, digits)
  }
  # Array is an ensemble; recursively row binding each matrix's dispersions
  else if(class(array) == "list"){
    purrr::map_dfr(array, .dispersion_matrix, pairs, sort_norms, singular, norm_pow, digits)
  }
}