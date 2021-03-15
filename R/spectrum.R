
# This script includes functions that help extract and plot eigenvalues of matrices.

#=================================================================================#
#                              SPECTRUM FUNCTIONS
#=================================================================================#

#' @title Obtain the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be returned
#' @param components returns the array with resolved real and imaginary components if TRUE; otherwise returns complex-valued eigenvalues
#' @param sortByNorm sorts the eigenvalue spectrum by its norms, otherwise sorts by sign
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
spectrum <- function(array, components = TRUE, sortByNorm = NA, order = NA){
  digits <- 4 # Digits to round values to
  sortByNorm <- .parse_sortByNorm(sortByNorm, array) # Parse for default values
  # Array is a matrix; call function returning eigenvalues for singleton matrix
  if(class(array) == "matrix"){.spectrum_matrix(array, components, sortByNorm, order, digits)}
  # Array is an ensemble; recursively row binding each matrix's eigenvalues
  else if(class(array) == "list"){purrr::map_dfr(array, .spectrum_matrix, components, sortByNorm, order, digits)}
}

# Helper function returning tidied eigenvalue array for a matrix
.spectrum_matrix <- function(P, components, sortByNorm, order, digits = 4){
  # Get the sorted eigenvalue spectrum of the matrix
  eigenvalues <- eigen(P)$values # Compute the eigenvalues of P
  if(sortByNorm){eigenvalues <- .sortByNorm(eigenvalues)} # Order the eigenvalue spectrum by norm rather than sign
  ## Filter for orders and evaluate spectrum
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

.parse_sortByNorm <- function(sortByNorm, array){
  if(class(array) == "list"){P <- array[[1]]} else{P <- array} # Parse array type to sample a matrix if ensemble
  # If the sortByNorm argument is uninitialized, infer optimal case. Optimally TRUE when eigenvalues are complex.
  if(is.na(sortByNorm)){sortByNorm <- ifelse(.isHermitian(P), F, T)} # Eigenvalues are real when the matrix is symmetric
  else{sortByNorm <- sortByNorm}
  sortByNorm # Return parsed value
}

#=================================================================================#
#                              ORDER SORTING SCHEMES
#=================================================================================#

# Sort an array of numbers by their norm (written for eigenvalue sorting)
.sortByNorm <- function(eigenvalues){(data.frame(eigenvalue = eigenvalues, norm = abs(eigenvalues)) %>% arrange(desc(norm)))$eigenvalue}

# Resorts the norm of an eigenvalue based on a particular metric of order; default is regular norm.
# .resort_spectrum <- function(array_spectrum, scheme = "norm"){
#   NA
  # Other schemes could be by sign; as in original "bug", we could always prioritize eigenvalues with positive sign
  # Could be called scheme = "sign"
  # array_spectrum
# }

#=================================================================================#
#                         SPECTRUM VISUALIZATION FUNCTIONS
#=================================================================================#

#' @title Visualize a plot of the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a scatterplot of the eigenvalues of a random matrix or ensemble.
#'
#' @inheritParams spectrum
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
  if(class(array) == "list" || class(array) == "matrix"){array_spectrum <- spectrum(array, ...)}
  else{array_spectrum <- array} # Else, the array is a precomputed spectrum (avoid computational waste for multiple visualizations)
  # Infer plot title string from which type of array (matrix/ensemble)
  title_str <- .plot_title(class(array), prefix = "Spectrum", mat_str)
  # Plot parameters
  order <- array_spectrum[["Order"]]
  # Plot
  array_spectrum %>%
    ggplot() +
    geom_point(mapping = aes(x = Re, y = Im, color = order), alpha = 0.75) +
    scale_color_continuous(type = "viridis") +
    labs(x = "Re", y = "Im", title = paste(title_str,sep = "")) +
    coord_fixed()
}

#' @title Visualize a plot of the eigenvalue distribution of a matrix or ensemble of matrices.
#'
#' @description Returns a histogram of the eigenvalues of a random matrix or ensemble.
#'
#' @inheritParams spectrum
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#' @param component a string specifying a specific component of the spectrum to display; either "Re" or "Im". Defaults to both
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
spectrum.histogram <- function(array, ..., component = NA, bins = 100, mat_str = ""){
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){array_spectrum <- spectrum(array, ...)}
  else{array_spectrum <- array} # Else, the array is a precomputed spectrum (avoid computational waste for multiple visualizations)
  # Infer plot title string from which type of array (matrix/ensemble)
  title_str <- .plot_title(class(array), prefix = "Spectrum", mat_str)
  # Plot parameters
  color0 <- "mediumpurple3"
  if(class(component) == "logical"){component <- c("Re", "Im")} # Set default to both components
  # Plot lambda function
  component_plot <- function(component){
    # Plot parameters
    component_str <- paste(" (",component,")",collapse = "")
    # Plot
    array_spectrum %>%
      ggplot() +
      geom_histogram(mapping = aes_string(x = component), fill = color0) +
      labs(x = component, y = "Frequency", title = paste(title_str, component_str, sep = ""))
  }
  # Get list of plots
  plots <- purrr::map(component, component_plot)
  # If we have both components and patchwork is loaded, attach plots to each other
  if(length(plots) == 2){plots[[1]] / plots[[2]]} else if(length(plots) == 1){plots[[1]]}
  # Return the list of plots
  else{plots}
}

# Helper function returning appoporiate string for matrix/ensemble given a matrix type string and class of input array
.plot_title <- function(array_class, prefix, mat_str){
  if(mat_str != ""){pre_space <- " "} else{pre_space <- ""} # Format without given name
  # Infer plot title string from which type of array (matrix/ensemble)
  if(array_class == "matrix"){plot_str <- paste(pre_space, mat_str," Matrix", sep = "", collapse = "")}
  else if(array_class == "list"){plot_str <- paste(pre_space, mat_str," Matrix Ensemble", sep = "", collapse = "")}
  else{plot_str <- paste(pre_space, mat_str," Matrix Ensemble", sep = "", collapse = "")}
  paste(prefix," of a",plot_str, sep = "")
}
