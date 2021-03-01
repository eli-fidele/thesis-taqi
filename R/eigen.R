
# This script includes functions that help extract and plot eigenvalues of matrices.

#=================================================================================#
#                              EIGENVALUE DISPERSION
#=================================================================================#

#' @title Obtain the eigenvalue spacings of a matrix or ensemble of matrices.
#'
#' @description Returns a vector of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalue spacings are to be returned
#' @param norm use the norm metric for eigenvalue spacing; otherwise returns absolute difference metric
#' @param components returns the array with resolved real and imaginary components; otherwise returns complex-valued vectors/distances
#' @param digits number of digits to round up values to
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
dispersion <- function(array, norm = T, components = T, digits = 3){
  # Array is a matrix; call function returning dispersion for singleton matrix
  if(class(array) == "matrix"){.dispersion_matrix(array, norm, components, digits)}
  # Array is an ensemble; recursively row binding each matrix's dispersions
  else if(class(array) == "list"){
    pairs <- .unique_pairs(nrow(array[[1]])) # Compute pairs to avoid computational waste and pass as argument
    purrr::map_dfr(.x = array, .f = .dispersion_matrix, norm, components, digits, pairs)
  }
}

# Find the eigenvalue dispersions for a given matrix
.dispersion_matrix <- function(P, norm = T, components = T, digits = 3, pairs = NA){
  eigenvalues <- eigen(P)$values # Get the eigenvalues of a matrix
  N <- nrow(P) # Get matrix dimension
  # If uninitialized for the ensemble, enumerate unique pairs of N eigenvalues
  if(class(pairs) == "logical"){idx_pairs <- .unique_pairs(N)} else{idx_pairs <- pairs} # Otherwise, read in pre-computed values
  # User is requesting a norm function rather than a setting of Euclidean norm
  if(class(norm) != "logical"){
    norm_fn <- function(x){(abs(x))^norm} # Beta norm
    purrr::map2_dfr(idx_pairs[,1], idx_pairs[,2], .resolve_dispersion, eigenvalues, norm_fn, components, digits)
  } 
  # Otherwise, just use the absolute value norm
  else{purrr::map2_dfr(idx_pairs[,1], idx_pairs[,2], .resolve_dispersion, eigenvalues, norm, components, digits)}
}

# Read and parse a dispersion observation between eigenvalue i and j.
.resolve_dispersion <- function(i, j, eigenvalues, norm, components, digits){
  # Compute the difference
  difference <- eigenvalues[i] - eigenvalues[j]
  # Resolve parameters of desired dispersion metric
  if(class(norm) == "function"){disp <- data.frame(Dispersion = norm(difference))}
  else if(norm){disp <- data.frame(Dispersion = abs(difference))}
  else{
    if(components){disp <- data.frame(Disp_Re = Re(difference), Disp_Im = Im(difference))}
    else{disp <- data.frame(Dispersion = difference)}
  }
  disp <- round(disp, digits) # Round digits
  cbind(disp, data.frame(OrderDiff = as.double(i - j))) # Add difference of order metric
}

# Enumerate the unique pairs given N items
.unique_pairs <- function(N){
  is <- do.call("c",map(1:N, function(i){rep(i,N)}))
  js <- rep(1:N, N)
  do.call("rbind",purrr::map2(is, js, .f = function(i, j){if(i > j){c(i = i, j = j)}}))
}
#=================================================================================#
#                         DISPERSION VISUALIZATION FUNCTIONS
#=================================================================================#

#' @title Visualize a plot of the eigenvalue difference spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a histogram of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues spacings are to be plotted
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
dispersion.histogram <- function(array, ..., bins = 100){
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, ...)}
  else{disps_df <- array} # Otherwise, the array is a precomputed dispersion dataframe
  num_entries <- nrow(disps_df) # Get number of entries
  # Plot parameters
  color0 <- "darkorchid4"
  # Return plot
  ggplot(data = disps_df, mapping = aes(x = Dispersion, y = stat(count / num_entries))) +
    geom_histogram(bins = bins) +
    labs(title = "Distribution of Eigenvalue Spacings", y = "Probability")
}

#' @title Visualize a plot of the eigenvalue difference spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a scatterplot of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues spacings are to be plotted
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
dispersion.scatterplot <- function(array, ...){
  # Process dispersion of the matrix/ensemble; if array is a dispersion data frame, copy it.
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, ...)}
  else{disps_df <- array} # Otherwise, the array is a precomputed dispersion dataframe
  # Plot parameters
  color0 <- "darkorchid4"
  # Get variances by level
  disps_df %>%
    ggplot(mapping = aes(x = Dispersion, y = OrderDiff, color = OrderDiff)) +
    geom_point() +
    scale_color_continuous(type = "viridis") +
    labs(title = "Distribution of Eigenvalue Spacings by Ranking Difference Class", y = "Ranking Difference")
}

.dispersion.varplot <- function(array, ...){
  # Process dispersion of the matrix/ensemble; if array is a dispersion data frame, copy it.
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, ...)}
  else{disps_df <- array}
  # Plot parameters
  color0 <- "darkorchid4"
  # Get variances by level
  disps_df %>%
    group_by(OrderDiff) %>%
    summarize(Var_Disp = var(Dispersion), size = n()) %>%
    ggplot(mapping = aes(x = OrderDiff, y = Var_Disp, color = Var_Disp, size = size)) +
    geom_point() +
    scale_color_continuous(type = "viridis") +
    #scale_size_manual(values = c(1,2)) +
    labs(title = "Variance of Eigenvalue Spacings by Ranking Difference Class", x = "Ranking Difference", y = "Variance")
}

#=================================================================================#
#                              SPECTRUM FUNCTIONS
#=================================================================================#

#' @title Obtain the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be returned
#' @param components returns the array with resolved real and imaginary components; otherwise returns complex-valued eigenvalues
#' @param largest returns the largest eigenvalues of the matrix (ensemble)
#' @param smallest returns the smallest eigenvalues of the matrix (ensemble)
#' @param digits number of digits to round up values to
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
spectrum <- function(array, components = T, largest = F, smallest = F, digits = 3){
  # Array is a matrix; call function returning eigenvalues for singleton matrix
  if(class(array) == "matrix"){.spectrum_matrix(array, components, largest, smallest, digits)}
  # Array is an ensemble; recursively row binding each matrix's eigenvalues
  else if(class(array) == "list"){purrr::map_dfr(array, .spectrum_matrix, components, largest, smallest, digits)}
}

# Helper function returning tidied eigenvalue array for a matrix
.spectrum_matrix <- function(P, components = T, largest = F, smallest = F, digits = 3){
  eigenvalues <- eigen(P)$values # Get eigenvalues of matrix P
  # Get largest eigenvalue
  if(largest){.resolve_eigenvalue(order = 1, eigenvalues, components)}
  # Get smallest eigenvalue
  else if(smallest){.resolve_eigenvalue(order = nrow(P), eigenvalues, components)}
  # Get all the eigenvalues
  else{purrr::map_dfr(1:nrow(P), .resolve_eigenvalue, eigenvalues, components, digits)}
}

# Read and parse an eigenvalue from an eigen(P)$value array
.resolve_eigenvalue <- function(order, eigenvalues, components, digits){
  # Read from eigen(P)$values
  eigenvalue <- eigenvalues[order]
  # If components are requested, resolve parts into seperate columns
  if(components){
    data.frame(Re = round(Re(eigenvalue),digits), Im = round(Im(eigenvalue),digits),
               Norm = round(abs(eigenvalue), digits), Order = order)
  }
  else{data.frame(Eigenvalue = round(eigenvalue, digits),
                  Norm = round(abs(eigenvalue), digits), Order = order)
  }
}

#=================================================================================#
#                         SPECTRUM VISUALIZATION FUNCTIONS
#=================================================================================#

#' @title Visualize a plot of the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a scatterplot of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be plotted
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#' @param mat_str (optional) a string argument of the class of the matrix to label the plot title.
#'
#' @return A ggplot object containing a scatterplot of the matrix/matrix ensemble's spectrum.
#' @examples
#' # Eigenvalue spectrum plot of a matrix
#' P <- RM_norm(N = 5)
#' #spectrum.scatterplot(P)
#'
#' # Labelled spectrum plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #spectrum.scatterplot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spectra plot of an ensemble of normal matrices
#' ensemble <- RME_norm(N = 3, size = 10)
#' #spectrum.scatterplot(ensemble)
#'
spectrum.scatterplot <- function(array, ..., mat_str = ""){
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){eigen_spectrum <- spectrum(array, ...)}
  else{eigen_spectrum <- array}
  # Infer plot title string from which type of array (matrix/ensemble)
  plot_str <- .plot_title(class(array), mat_str)
  # Plot parameters
  order <- eigen_spectrum[["Order"]]
  # Plot
  eigen_spectrum %>%
    ggplot() +
    geom_point(mapping = aes(x = Re, y = Im, color = order), alpha = 0.75) +
    scale_color_continuous(type = "viridis") +
    labs(x = "Re", y = "Im", title = paste("Spectrum of a",plot_str,sep = "")) +
    coord_fixed()
}

#' @title Visualize a plot of the eigenvalue distribution of a matrix or ensemble of matrices.
#'
#' @description Returns a histogram of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be plotted
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#' @param imaginary if TRUE, returns the distribution of imaginary components
#' @param bins number of bins of the histogram
#' @param mat_str (optional) a string argument of the class of the matrix to label the plot title.
#'
#' @return A ggplot object containing a histogram of the matrix/matrix ensemble's spectrum.
#' @examples
#' # Eigenvalue spectrum plot of a matrix
#' P <- RM_norm(N = 5)
#' #spectrum.histogram(P)
#'
#' # Labelled spectrum plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #spectrum.histogram(Q, mat_str = "Beta")
#'
#' # Eigenvalue spectra plot of an ensemble of normal matrices
#' ensemble <- RME_norm(N = 3, size = 10)
#' #spectrum.histogram(ensemble)
#'
spectrum.histogram <- function(array, ..., imaginary = F, bins = 100, mat_str = ""){
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){eigen_spectrum <- spectrum(array, ...)}
  else{eigen_spectrum <- array}
  # Infer plot title string from which type of array (matrix/ensemble)
  plot_str <- .plot_title(class(array), mat_str)
  # Plot parameters
  color0 <- "darkorchid4"
  if(imaginary){component <- "Im"} else{component <- "Re"}
  # Plot
  eigen_spectrum %>%
    ggplot() +
    geom_histogram(mapping = aes_string(x = component), fill = color0) +
    labs(x = component, y = "Frequency", title = paste("Spectrum of a",plot_str,sep = ""))
}

# Helper function returning appoporiate string for matrix/ensemble given a matrix type string and class of input array
.plot_title <- function(array_class, mat_str){
  if(mat_str != ""){pre_space <- " "} else{pre_space <- ""} # Format without given name
  # Infer plot title string from which type of array (matrix/ensemble)
  if(array_class == "matrix"){plot_str <- paste(pre_space, mat_str," Matrix", sep = "", collapse = "")}
  else if(array_class == "list"){plot_str <- paste(pre_space, mat_str," Matrix Ensemble", sep = "", collapse = "")}
  plot_str
}