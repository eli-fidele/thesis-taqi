
# This script includes functions that help extract and plot eigenvalues of matrices.

#=================================================================================#
#                              SPECTRUM FUNCTIONS
#=================================================================================#

#' @title Obtain the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be returned
#' @param order get eigenvalues with that given order (norm ranking); order 1 represents largest, order N represents smallest (where N is the number of eigenvalues).
#'   If uninitialized, returns the entire spectrum.
#' @param components returns the array with resolved real and imaginary components if TRUE; otherwise returns complex-valued eigenvalues
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
spectrum <- function(array, order = NA, components = T, digits = 3){
  # Array is a matrix; call function returning eigenvalues for singleton matrix
  if(class(array) == "matrix"){.spectrum_matrix(array, order, components, digits)}
  # Array is an ensemble; recursively row binding each matrix's eigenvalues
  else if(class(array) == "list"){purrr::map_dfr(array, .spectrum_matrix, order, components, digits)}
}

# Helper function returning tidied eigenvalue array for a matrix
.spectrum_matrix <- function(P, order = NA, components = T, digits = 3){
  eigenvalues <- eigen(P)$values # Get eigenvalues of matrix P
  if(class(order) == "logical"){order <- 1:nrow(P)} # If uninitialized, get eigenvalues of all orders
  else{order <- c(order)} # Otherwise, concatenate so single inputs become vectors
  purrr::map_dfr(order, .resolve_eigenvalue, eigenvalues, components, digits) # Get the eigenvalues
}

# Read and parse an eigenvalue from an eigen(P)$value array
.resolve_eigenvalue <- function(order, eigenvalues, components, digits){
  eigenvalue <- eigenvalues[order] # Read from eigen(P)$values
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
#' @inheritParams spectrum
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
