
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
spectrum <- function(array, order = NA, components = T, digits = 3){ #### Add sortBy = "norm" or "sign" feature????
  # Array is a matrix; call function returning eigenvalues for singleton matrix
  if(class(array) == "matrix"){.spectrum_matrix(array, order, components, digits)}
  # Array is an ensemble; recursively row binding each matrix's eigenvalues
  else if(class(array) == "list"){purrr::map_dfr(array, .spectrum_matrix, order, components, digits)}
}

# Helper function returning tidied eigenvalue array for a matrix
.spectrum_matrix <- function(P, order = NA, components = T, digits = 3){
  eigenvalues <- .sort_norm(eigen(P)$values) # Get eigenvalues of matrix P, sorted by size
  if(class(order) == "logical"){order <- 1:nrow(P)} # If uninitialized, get eigenvalues of all orders
  else{order <- c(order)} # Otherwise, concatenate so single inputs become vectors
  purrr::map_dfr(order, .resolve_eigenvalue, eigenvalues, components, digits) # Get the eigenvalues
}

# Parses a matrix spectrum array for the eigenvalue at a given order as cplx type (for arithmetic)
.read_eigenvalue <- function(order, mat_spectrum){
  if(ncol(mat_spectrum) == 3){mat_spectrum[order, 1]} # If the components are resolved, return value in the first (Eigenvalue) column
  else{ # Components are resolved; get components and make it a complex number for arithmetic prep
    evalue <- complex(real = mat_spectrum[order, 1], imaginary = mat_spectrum[order, 2])
    if(Im(evalue) != 0){evalue} else{as.numeric(evalue)} # If it is real, coerce it into a numeric to remove +0i 
    }
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

# Resorts the norm of an eigenvalue based on a particular metric of order; default is regular norm.
resort_spectrum <- function(array_spectrum, scheme = "norm"){
  NA
  # Other schemes could be by sign; as in original "bug", we could always prioritize eigenvalues with positive sign
  # Could be called scheme = "sign"
  array_spectrum
}

# Sort an array of numbers by their norm (written for eigenvalue sorting)
.sort_norm <- function(eigenvalues){
  (data.frame(eigenvalue = eigenvalues, norm = abs(eigenvalues)) %>% arrange(desc(norm)))$eigenvalue
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
  plots <- map(component, component_plot)
  # If we have both components and patchwork is loaded, attach plots to each other
  if(require(patchwork) && length(plots) == 2){plots[[1]] / plots[[2]]} else if(length(plots) == 1){plots[[1]]}
  # Return the list of plots
  else{plot_list}
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
