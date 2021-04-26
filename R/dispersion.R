
# dispersion : A script file including functions to help compute the dispersions
# of random matrices and random matrix ensembles.

#=================================================================================#
#                                   DISPERSION
#=================================================================================#

#' @title Obtain the eigenvalue spacings of a matrix or ensemble of matrices.
#' @description Returns a vector of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalue spacings are to be returned
#' @param pairs a string argument representing the pairing scheme to use
#' @param norm_order sorts the eigenvalue spectrum by its norms when TRUE; otherwise, sorts eigenvalue by sign
#' @param singular get the singular values of the matrix (i.e. square root of the eigenvalues of the matrix times its transpose)
#' @param pow_norm power to raise norm to - defaults to 1 (the standard absolute value); otherwise raises norm to the power of argument (beta norm)
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#'
#' @examples
#' # Eigenvalue dispersion of a normal matrix using the lower pair scheme
#' P <- RM_norm(N = 5)
#' disp_P <- dispersion(P, pairs = "lower")
#'
#' # Eigenvalue dispersion of a stochastic matrix (using the consecutive pair scheme)
#' Q <- RM_stoch(N = 5)
#' disp_Q <- dispersion(Q)
#'
#' # Eigenvalue dispersion of an normal matrix ensemble, ordering by sign instead of norm.
#' ens <- RME_beta(N = 10, beta = 2, size = 10)
#' disp_ens <- dispersion(ens, norm_order = FALSE)
#'
dispersion <- function(array, pairs = NA, norm_order = TRUE, singular = FALSE, pow_norm = 1){
  # Digits to round values to
  digits <- 4
  # Get the type of array
  array_class <- .arrayClass(array)
  # Parse input and generate pair scheme (default NA), passing on array for dimension
  pairs <- .parsePairs(pairs, array, array_class)
  # Array is an ensemble; recursively row binding each matrix's dispersions
  if(array_class == "ensemble"){
    disp <- purrr::map_dfr(array, .dispersion_matrix, pairs, norm_order, singular, pow_norm, digits)
  }
  # Array is a matrix; call function returning dispersion for singleton matrix
  else if(array_class == "matrix"){
    disp <- .dispersion_matrix(array, pairs, norm_order, singular, pow_norm, digits)
  }
  # Resolve column types; i.e. coerce real-valued eigenvalues to a numeric type if possible
  disp <- .resolveNumType(disp)
  # Return the dispersion
  disp
}

#=================================================================================#
# Find the eigenvalue dispersions for a given matrix
.dispersion_matrix <- function(P, pairs, norm_order, singular, pow_norm, digits = 4){
  # Get the ordered spectrum of the matrix
  eigenvalues <- spectrum(P, norm_order = norm_order, singular = singular)
  # Generate norm function to pass along as argument (Euclidean or Beta norm)
  norm_fn <- function(x){ (abs(x))^pow_norm }
  # Compute the dispersion
  disp <- purrr::map2_dfr(pairs[["i"]], pairs[["j"]], .resolve_dispersion, eigenvalues, norm_fn, digits)
  # Return the dispersion
  disp
}

#=================================================================================#
# Read and parse a dispersion observation between eigenvalue i and j.
.resolve_dispersion <- function(i, j, eigenvalues, norm_fn, digits){
  ## Copmute dispersion metrics
  disp <- data.frame(i = i, j = j) # Initialize dispersion dataframe by adding order of eigenvalues compared
  disp$eig_i <- .read_eigenvalue(i, eigenvalues); disp$eig_j <- .read_eigenvalue(j, eigenvalues) # Add the eigenvalues
  disp$id_diff <- disp$eig_j - disp$eig_i # Get the identity difference dispersion metric
  ## Compute norm metrics
  disp$id_diff_norm <- norm_fn(disp$id_diff) # Take the norm of the difference
  disp$abs_diff <- norm_fn(disp$eig_j) - norm_fn(disp$eig_i) # Compute the difference of absolutes w.r.t. norm function (Euclidean or beta)
  ## Prepare for return
  disp <- round(disp, digits) # Round digits
  disp$diff_ij <- disp$i - disp$j
  disp # Return resolved dispersion observation
}

#=================================================================================#
#                           DISPERSION: HELPER FUNCTIONS
#=================================================================================#

# Parses a matrix spectrum array for the eigenvalue at a given order as cplx type (for arithmetic)
.read_eigenvalue <- function(order, mat_spectrum){
  # If the components are not resolved, return value in the first (Eigenvalue) column
  if(ncol(mat_spectrum) == 3){ mat_spectrum[order, "Eigenvalue"] }
  # Components are resolved; get components and make it a complex number for arithmetic prep
  else{ complex(real = mat_spectrum[order, "Re"], imaginary = mat_spectrum[order, "Im"]) }
}

#=================================================================================#
# Resolves the numerical types of the eigenvalue columns of the dispersion dataframe
.resolveNumType <- function(disp){
  # See if the eigenvalues are complex
  is_complex <- FALSE %in% c(Im(disp$eig_i) == 0, Im(disp$eig_j) == 0)
  # If it is real, coerce it into a numeric, removing the +0i
  if(!is_complex){
    disp$eig_i <- as.numeric(disp$eig_i)
    disp$eig_j <- as.numeric(disp$eig_j)
    disp$id_diff <- as.numeric(disp$id_diff)
  }
  disp # Return the resolved dispersion
}

#=================================================================================#
# Parse a string argument for which pairing scheme to utilize
.parsePairs <- function(pairs, array, array_class){
  valid_schemes <- c("largest", "lower", "upper", "consecutive", "all") # Valid schemes for printing if user is unaware of options
  # Obtain the matrix by inferring array type; if ensemble take first matrix
  if(array_class == "ensemble"){
    P <- array[[1]]
  } else if(array_class == "matrix"){
    P <- array
  }
  if(class(pairs) == "logical"){pairs <- "consecutive"} # Set default value to be the consecutive pair scheme
  # Stop function call if the argument is invalid
  if(!(pairs %in% valid_schemes)){
    scheme_list <- paste(valid_schemes, collapse = ", ")
    stop(paste("Invalid pair scheme. Try one of the following: ", scheme_list, ".", ""))
  }
  # Obtain the dimension of the matrix
  N <- nrow(P)
  # Parse the pair string and evaluate the pair scheme
  if(pairs == "largest"){pair_scheme <- data.frame(i = 2, j = 1)}
  else if(pairs == "consecutive"){pair_scheme <- .consecutive_pairs(N)}
  else if(pairs == "lower"){pair_scheme <- .unique_pairs_lower(N)}
  else if(pairs == "upper"){pair_scheme <- .unique_pairs_upper(N)}
  else if(pairs == "all"){pair_scheme <- .all_pairs(N)}
  pair_scheme # Return pair scheme
}

#=================================================================================#
#                                PAIRING SCHEMES
#=================================================================================#

# The trivial pairing scheme:
# Enumerate all possible pairs.
.all_pairs <- function(N){
  purrr::map_dfr(1:N, function(i, N){data.frame(i = rep(i, N), j = 1:N)}, N)
}

#=================================================================================#
# The consecutive pairing scheme:
# Enumerate all possible consecutive/neighboring pairs. No linear combiantions of dispersions exist.
.consecutive_pairs <- function(N){
  purrr::map_dfr(2:N, function(i){data.frame(i = i, j = as.integer(i - 1))})
}

#=================================================================================#
# The lower-triangular pairing scheme:
# Enumerate the pair combinations given N items with i > j.
.unique_pairs_lower <- function(N){
  is <- do.call("c", purrr::map(1:N, function(i){rep(i,N)}))
  js <- rep(1:N, N)
  pairs <- do.call("rbind",purrr::map2(is, js, .f = function(i, j){if(i > j){c(i = i, j = j)}}))
  data.frame(pairs)
}

#=================================================================================#
# The upper-triangular pairing scheme:
# Enumerate the pair combinations given N items with i < j.
.unique_pairs_upper <- function(N){
  is <- do.call("c", purrr::map(1:N, function(i){rep(i,N)}))
  js <- rep(1:N, N)
  pairs <- do.call("rbind",purrr::map2(is, js, .f = function(i, j){if(i < j){c(i = i, j = j)}}))
  data.frame(pairs)
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
  pairs <- .parsePairs(pairs, array)
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

