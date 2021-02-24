
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
spectrum <- function(array, largest = F, smallest = F){
  # Infer type of array (matrix or ensemble) then parse accordingly.
  is_ensemble <- (class(array) == "list")
  # One type of array is inferred, obtain the eigenvalue array
  if(!is_ensemble){
    P <- array 
    eigen_array <- data.frame(eigen(P)$values) # Get eigenvalues of matrix
    spectrum_row <- function(i, tbl){c(round(Re(tbl[i,]), 5), round(Im(tbl[i,]), 5), abs(tbl[i,]), i)}
    eigenvalues <- do.call("rbind", lapply(X = 1:nrow(P), FUN = spectrum_row, tbl = eigen_array)) # Extract eigenvalues
    eigenvalues <- data.frame(eigenvalues) # Array of eigenvalues
    colnames(eigenvalues) <- c("Re", "Im", "Norm", "Order") # Rename columns
  }
  # Otherwise, recursively obtain the ensemble's spectrum by row binding each matrix's returned spectrum
  else{
    ensemble <- array
    eigenvalues <- .ensemble_spectrum(ensemble)
  }
  # Once the eigenvalue array is obtained, filter for wanted statistics
  if(smallest){eigenvalues <- eigenvalues[which(eigenvalues$Order == which.min(eigenvalues$Order)),]}
  if(largest){eigenvalues <- eigenvalues[which(eigenvalues$Order == 1),]}
  #if(class(order) != NULL){eigenvalues <- eigenvalues %>% filter(Order %in% c(order))}
  # Return spectrum of eigenvalues
  eigenvalues
}

# Returns largest eigenvalues of the matrix
.LRGST <- function(spectrum){}
.SMLST <- function(spectrum){spectrum[which(spectrum$Index == 1)]}

# Helper function for spectrum, returns a tidied dataframe of the eigenvalues of a matrix ensemble input
.ensemble_spectrum <- function(ensemble){
  do.call("rbind",lapply(X = 1:length(ensemble), FUN = function(i){spectrum(ensemble[[i]])})) # Return the spectra for this ensemble
}

# Read in the eigenvalue in the Kth row from a eigenvalue array and return a numerical
.read_eigenvalue <- function(spectrum, K){
  # if im is 0 maybe return real val
  # filter for index (largest == 1)
  complex(real = spectrum[K,1], imaginary = spectrum[K,2])
  }

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
spectrum.plot <- function(array, mat_str = ""){
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){eigen_spectrum <- spectrum(array)}
  else{eigen_spectrum <- array}
  # Infer plot title string from which type of array (matrix/ensemble)
  is_mat <- class(array) == "matrix"
  if(is_mat){class_str <- paste("Matrix", sep = "")}
  else{class_str <- paste("Matrix Ensemble", sep = "")}
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
    geom_point(mapping = aes(x = Re, y = Im, color = Order), alpha = 0.75) +
    scale_color_continuous(type = "viridis") +
    #theme(legend.position = "none") +
    labs(x = "Re", y = "Im", title = paste("Spectrum of a",mat_str,class_str,sep = " ")) #+
    #xlim(x_range) +
    #ylim(-r,r) +
    #coord_fixed(ratio = 1)
}

