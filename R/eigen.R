
# This script includes functions that help extract and plot eigenvalues of matrices.


#=================================================================================#
#                              SPECTRUM FUNCTIONS
#=================================================================================#

#' @title Obtain the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a tidied dataframe of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be returned
#'
#' @return A tidy dataframe with the real & imaginary components of the eigenvalues and their norms along with a unique index.
#' @examples
#' # Eigenvalue spectrum of a matrix
#' P <- RM_norm(N = 5)
#' spectrum_P <- spectrum(P)
#'
#' Q <- matrix(runif(2^2), ncol = 2)
#' spectrum_Q <- spectrum(Q)
#'
#' # Eigenvalue spectra of ensemble matrices
#' ensemble <- RME("norm", args = c(N = 3), ensemble_size = 10)
#' ensemble_spectrum <- spectrum(ensemble)
#'
spectrum <- function(array){
  # See if we have a ensemble of matrices or a single matrix
  is_ensemble <- (class(array) == "list")
  # If input is a matrix, proceed to get its spectrum
  if(!is_ensemble){
    P <- array # Rename array
    N <- nrow(P) # Obtain matrix dimension
    eigen_array <- data.frame(eigen(P)$values) # Get eigenvalues
    spectrum_row <- function(i, tbl){c(round(Re(tbl[i,]), 5), round(Im(tbl[i,]), 5), abs(tbl[i,]), i)}
    eigenvalues <- do.call("rbind", lapply(X = 1:N, FUN = spectrum_row, tbl = eigen_array)) # Extract eigenvalues
    colnames(eigenvalues) <- c("Re", "Im", "Norm", "Index") # Rename columns
    return(data.frame(eigenvalues))
  }
  # Otherwise, recursively obtain the ensemble's spectrum by row binding each matrix's returned spectrum
  else{
    ensemble <- array # Rename array
    return(.ensemble_spectrum(ensemble))
  }
}

# Helper function for spectrum, returns a tidied dataframe of the eigenvalues of a matrix ensemble input
.ensemble_spectrum <- function(ensemble){
  do.call("rbind",lapply(X = 1:length(ensemble), FUN = function(i){spectrum(ensemble[[i]])})) # Return the spectra for this ensemble
}

# Read in the eigenvalue in the Kth row from a eigenvalue array and return a numerical
.read_eigenvalue <- function(spectrum, K){complex(real = spectrum[K,1], imaginary = spectrum[K,2])}

#=================================================================================#
#                         SPECTRUM VISUALIZATION FUNCTIONS
#=================================================================================#

#' @title Visualize a plot of the eigenvalue spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a scatterplot of the eigenvalues of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues are to be plotted
#' @param mat_str (optional) a string argument of the class of the matrix to label the plot title.
#'
#' @return A ggplot object containing a scatterplot of the matrix/matrix ensemble's spectrum.
#' @examples
#' # Eigenvalue spectrum of a matrix
#' P <- RM_norm(N = 5)
#' spectrum_plot(P)
#'
#' # Labelled spectrum plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' spectrum_plot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spectra of ensemble matrices
#' ensemble <- RME("norm", args = c(N = 3), ensemble_size = 10)
#' spectrum_plot(ensemble)
#'
spectrum_plot <- function(array, mat_str = ""){
  # See if we have a ensemble of matrices or a single matrix
  is_ensemble <- (class(array) != "matrix")
  # If not ensemble, directly process the spectrum of the matrix
  if(!is_ensemble){
    P <- array # Rename matrix
    eigen_spectrum <- spectrum(P)
    mat_str <- paste(mat_str, "Matrix", sep = " ")
    }
  else{ # Otherwise, process the ensemble first, and update the plot string.
    ensemble <- array
    eigen_spectrum <- .ensemble_spectrum(ensemble)
    mat_str <- paste(mat_str, "Matrix Ensemble", sep = " ")
    }
  # Plot parameters
  #r <- 1
  #x_window <- 0.5
  #x_range <- c(-(r + x_window), (r + x_window)) # Widen the width of the plot
  #circle <- data.frame(x0 = 0, y0 = 0, r = r)
  # Color plot parameters
  color0 <- "steelblue"
  color1 <- "deepskyblue3"
  # Plot
  ggplot2::ggplot(eigen_spectrum) +
    #geom_circle(mapping = aes(x0 = x0, y0 = y0, r = r), data = circle, color = color0) +
    geom_point(mapping = aes(x = Re, y = Im, color = Re), alpha = 0.75) +
    scale_color_continuous(type = "viridis") +
    theme(legend.position = "none") +
    labs(x = "Re", y = "Im", title = paste("Spectrum of a ",mat_str,sep = "")) #+
    #xlim(x_range) +
    #ylim(-r,r) +
    #coord_fixed(ratio = 1)
}

