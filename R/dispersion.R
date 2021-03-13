
# This script includes functions that help extract and plot eigenvalues of matrices.

#=================================================================================#
#                              EIGENVALUE DISPERSION
#=================================================================================#

#' @title Obtain the eigenvalue spacings of a matrix or ensemble of matrices.
#'
#' @description Returns a vector of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalue spacings are to be returned
#' @param pairs a string argument representing the pairing scheme to use 
#' @param sortByNorm sorts the eigenvalue spectrum by its norms when TRUE; otherwise, sorts eigenvalue by sign
#' @param norm_pow power to raise norm to - defaults to 1 (the standard absolute value); otherwise raises norm to the power of argument (beta norm)
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#' @examples
#'
#' # Eigenvalue dispersion of a normal matrix
#' P <- RM_norm(N = 5)
#' #disp_P <- dispersion(P)
#'
#' # Eigenvalue dispersion of a stochastic matrix
#' Q <- RM_stoch(N = 5)
#' #disp_Q <- dispersion(Q)
#'
#' # Eigenvalue dispersion of an ensemble
#' ensemble <- RME_norm(N = 3, size = 10)
#' #disp_ensemble <- dispersion(ensemble)
#'
#' # Alternatively, use the pipe
#' #disp_ensemble <- RME_norm(N = 3, size = 10) %>% dispersion()
#'
dispersion <- function(array, pairs = NA, sortByNorm = NA, norm_pow = 1){ #sortNorms? orderByNorms? pair_scheme?
  digits <- 4 # Digits to round values to
  ## Initialize pair scheme by parsing input ##
  pairs <- .parse_pairs(pairs, P) # Generate pairs by parsing the string input (or default NA), passing on array for dimension and array type inference
  # Array is a matrix; call function returning dispersion for singleton matrix
  if(class(array) == "matrix"){.dispersion_matrix(array, pairs, sortByNorm, norm_pow, digits)}
  # Array is an ensemble; recursively row binding each matrix's dispersions
  else if(class(array) == "list"){purrr::map_dfr(array, .dispersion_matrix, pairs, sortByNorm, norm_pow, digits)} 
}

# Find the eigenvalue dispersions for a given matrix
.dispersion_matrix <- function(P, pairs, sortByNorm, norm_pow, digits = 4){
  ## Get the sorted eigenvalue spectrum of the matrix ##
  eigenvalues <- spectrum(P, sortByNorm = sortByNorm)
  ## Get and evaluate the matrix dispersion ##
  N <- nrow(P) # Get matrix dimension
  norm_fn <- function(x){(abs(x))^norm_pow} # Generate norm function to pass along as argument (Euclidean or Beta norm)
  purrr::map2_dfr(pairs[,1], pairs[,2], .resolve_dispersion, eigenvalues, norm_fn, digits)
}

# Read and parse a dispersion observation between eigenvalue i and j.
.resolve_dispersion <- function(i, j, eigenvalues, norm_fn, digits){
  ## Copmute dispersion metrics
  disp <- data.frame(i = i, j = j) # Initialize dispersion dataframe by adding order of eigenvalues compared
  disp$eig_i <- .read_eigenvalue(i, eigenvalues); disp$eig_j <- .read_eigenvalue(j, eigenvalues) # Add the eigenvalues
  disp$id_diff <- disp$eig_j - disp$eig_i # Get the identity difference dispersion metric
  disp$id_diff_norm <- norm_fn(disp$id_diff) # Take the norm of the difference
  disp$abs_diff <- norm_fn(disp$eig_j) - norm_fn(disp$eig_i) # Compute the difference of absolutes w.r.t. norm function (Euclidean or beta)
  ## Prepare for return
  disp <- round(disp, digits) # Round digits
  disp$orderDiff_ji <- disp$i - disp$j
  disp # Return resolved dispersion observation
}

# Parses a matrix spectrum array for the eigenvalue at a given order as cplx type (for arithmetic)
.read_eigenvalue <- function(order, mat_spectrum){
  if(ncol(mat_spectrum) == 3){mat_spectrum[order, 1]} # If the components are resolved, return value in the first (Eigenvalue) column
  else{ # Components are resolved; get components and make it a complex number for arithmetic prep
    evalue <- complex(real = mat_spectrum[order, 1], imaginary = mat_spectrum[order, 2])
    if(Im(evalue) != 0){evalue} else{as.numeric(evalue)} # If it is real, coerce it into a numeric to remove +0i 
  }
}

#=================================================================================#
#                                 PAIR SCHEMA
#=================================================================================#

# Parse a string argument for which pairing scheme to utilize
.parse_pairs <- function(pairs, array){
  valid_schemes <- c("largest", "lower", "upper", "consecutive") # Valid schemes for printing if user is unaware of options
  # Obtain the matrix by inferring array type; if ensemble take first matrix
  if(class(array) == "list"){P <- array[[1]]} else{P <- array}
  if(class(pairs) == "logical"){pairs <- "consecutive"} # Set default value to be the consecutive pair scheme
  # Stop function call if the argument is invalid
  if(!(pairs %in% valid_schemes)){stop(paste("Invalid pair scheme. Try one of the following: ",paste(valid_schemes, collapse = ", "),".", sep = ""))} 
  # Obtain the dimension of the matrix
  N <- nrow(P)
  # Parse the pair string and evaluate the pair scheme
  if(pairs == "largest"){pair_scheme <- data.frame(i = 1, j = 2)}
  else if(pairs == "consecutive"){pair_scheme <- .consecutive_pairs(N)}
  else if(pairs == "lower"){pair_scheme <- .unique_pairs_lower(N)}
  else if(pairs == "upper"){pair_scheme <- .unique_pairs_upper(N)}
  pair_scheme # Return pair scheme
}

# The antisymmetric pair scheme (for assymetric dispersion metrics); essentially all the permutations  
.all_pairs <- function(N){
  purrr::map_dfr(1:N, function(i, N){data.frame(i = rep(i, N), j = 1:N)}, N)
}

# The consecutive-value scheme (Sufficient such that no linear combiantions of the diseprsion metric exists); one apart
.consecutive_pairs <- function(N){
  map_dfr(2:N, function(i){data.frame(i = i, j = i - 1)})
}

# The triangular pair schema (for symmetric dispersion metrics); essentially all the combinations
# (Enumerate the pair combinations given N items with i > j)
.unique_pairs_lower <- function(N){
  is <- do.call("c",map(1:N, function(i){rep(i,N)}))
  js <- rep(1:N, N)
  do.call("rbind",purrr::map2(is, js, .f = function(i, j){if(i > j){c(i = i, j = j)}}))
}
# The triangular pair schema (for symmetric dispersion metrics); essentially all the combinations
# (Enumerate the pair combinations given N items with i < j)
.unique_pairs_upper <- function(N){
  is <- do.call("c",map(1:N, function(i){rep(i,N)}))
  js <- rep(1:N, N)
  do.call("rbind",purrr::map2(is, js, .f = function(i, j){if(i < j){c(i = i, j = j)}}))
}
#=================================================================================#
#                         DISPERSION VISUALIZATION FUNCTIONS
#=================================================================================#

#' @title Visualize a plot of the eigenvalue difference spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a histogram of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#' @param metric a string denoting the eigenvalue dispersion metric (column of a dispersion object) to use
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#' @param bins (optional) a string argument of the class of the matrix to label the plot title.
#'
#' @return A ggplot object containing a histogram of the matrix/matrix ensemble's eigenvalue spacings.
#' @examples
#' # Eigenvalue spacings plot of a normal matrix
#' P <- RM_norm(N = 5)
#' #dispersion.histogram(P)
#'
#' # Eigenvalue spacings plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #dispersion.histogram(Q, mat_str = "Beta")
#'
#' # Eigenvalue spacings plot of an ensemble of normal matrices
#' # ensemble <- RME_norm(N = 3, size = 10)
#' # dispersion.histogram(ensemble)
#'
dispersion.histogram <- function(array, metric = NA, ..., bins = 100){
  valid_schemes <- c("id_diff","id_diff_norm","abs_diff") # Valid schemes for printing if user is unaware of options
  if(class(metric) == "logical"){stop("Please input a valid dispersion metric. Try one of the following: ",paste(valid_schemes, collapse = ", "),".", sep = "")}
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, ...)}
  else{disps_df <- array} # Otherwise, the array is a precomputed dispersion dataframe
  num_entries <- nrow(disps_df) # Get number of entries
  # Plot parameters
  color0 <- "darkorchid4"
  # Return plot
  ggplot(data = disps_df, mapping = aes_string(x = metric)) +
    geom_histogram(mapping = aes(y = stat(count / num_entries)), bins = bins) +
    labs(title = "Distribution of Eigenvalue Spacings", y = "Probability")
}

#' @title Visualize a plot of the eigenvalue difference spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a scatterplot of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#' @param metric a string denoting the eigenvalue dispersion metric (column of a dispersion object) to use
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#'
#' @return A ggplot object containing a scatterplot of the matrix/matrix ensemble's eigenvalue spacings.
#' @examples
#' # Eigenvalue spacings plot of a normal matrix
#' P <- RM_norm(N = 5)
#' #dispersion.scatterplot(P)
#'
#' # Eigenvalue spacings plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #dispersion.scatterplot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spacings plot of an ensemble of normal matrices
#' # ensemble <- RME_norm(N = 3, size = 10)
#' # dispersion.scatterplot(ensemble)
#'
dispersion.scatterplot <- function(array, metric = "id_diff_norm", ...){
  # Process dispersion of the matrix/ensemble; if array is a dispersion data frame, copy it.
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, ...)}
  else{disps_df <- array} # Otherwise, the array is a precomputed dispersion dataframe
  # Plot parameters
  color0 <- "darkorchid4"
  # Get variances by level
  disps_df %>%
    ggplot(mapping = aes(y = j, color = j)) + # If we have consecutive pairs, more useful. Add logic to parse for this feature
    #ggplot(mapping = aes(y = orderDiff_ji, color = orderDiff_ji)) +
    geom_point(mapping = aes_string(x = metric)) +
    scale_color_continuous(type = "viridis") +
    labs(title = "Distribution of Eigenvalue Spacings by Ranking Difference Class", y = "Ranking Difference")
}

#' @title Visualize a plot of the variances of the eigenvalue dispersions by order difference level given a matrix or an ensemble.
#'
#' @description Returns a variance scatterplot of classes/levels of eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#' @param metric a string denoting the eigenvalue dispersion metric (column of a dispersion object) to use
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#'
#' @return A ggplot object containing a scatterplot of the matrix/matrix ensemble's eigenvalue spacings.
#' @examples
#' # Eigenvalue spacings plot of a normal matrix
#' P <- RM_norm(N = 5)
#' #.dispersion.varplot(P)
#'
#' # Eigenvalue spacings plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #.dispersion.varplot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spacings plot of an ensemble of normal matrices
#' # ensemble <- RME_norm(N = 3, size = 10)
#' # .dispersion.varplot(ensemble)
#'
.dispersion.varplot <- function(array, metric, ...){
  # Process dispersion of the matrix/ensemble; if array is a dispersion data frame, copy it.
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, ...)}
  else{disps_df <- array}
  # Plot parameters
  color0 <- "darkorchid4"
  # Get variances by level
  disps_df %>%
    group_by(orderDiff_ji) %>%
    summarize(Var_Disp = var(Dispersion), size = n()) %>%
    ggplot(mapping = aes(x = orderDiff_ji, y = Var_Disp, color = Var_Disp, size = size)) +
    geom_point() +
    scale_color_continuous(type = "viridis") +
    #scale_size_manual(values = c(1,2)) +
    labs(title = "Variance of Eigenvalue Spacings by Ranking Difference Class", x = "Ranking Difference", y = "Variance")
}
