
# This script includes functions that help extract and plot eigenvalues of matrices.


#=================================================================================#
#                              SPECTRUM FUNCTIONS
#=================================================================================#

# Returns a tidied dataframe of the eigenvalues of a matrix ensemble
ensemble_spectrum <- function(ensemble){
  do.call("rbind",lapply(X = 1:length(ensemble), FUN = function(i){spectrum(ensemble[[i]])})) # Return the spectra for this ensemble
}

# Returns a tidied dataframe of the eigenvalues of a random matrix
spectrum <- function(P){
  N <- nrow(P) # Obtain matrix dimension
  eigen_array <- data.frame(eigen(P)$values) # Get eigenvalues
  spectrum_row <- function(i, array){c(round(Re(array[i,]), 5), round(Im(array[i,]), 5), abs(array[i,]), i)}
  eigenvalues <- do.call("rbind", lapply(X = 1:N, FUN = spectrum_row, eigen_array))
  colnames(eigenvalues) <- c("Re", "Im", "Norm", "Index") # Rename columns
  data.frame(eigenvalues)
}

#=================================================================================#
#                         SPECTRUM VISUALIZATION FUNCTIONS
#=================================================================================#

# Plots the eigenvalues of a given matrix P
spectrum_plot <- function(P, mat_str = ""){
  # See if we have a ensemble of matrices or a single matrix
  not_ensemble <- (class(nrow(P) == ncol(P)) == "numeric")
  # If not ensemble, directly process the spectrum
  if(not_ensemble){
    eigen_spectrum <- spectrum(P)
    mat_str <- paste(mat_str, "Matrix", sep = " ")
    }
  else{ # Otherwise, process the ensemble first, and update the plot string.
    eigen_spectrum <- ensemble_spectrum(P)
    mat_str <- paste(mat_str, "Matrix Ensemble", sep = " ")
    }
  # Plot parameters
  r <- 1
  x_window <- 0.5
  x_range <- c(-(r + x_window), (r + x_window)) # Widen the width of the plot
  circle <- data.frame(x0 = 0, y0 = 0, r = r)
  # Color plot parameters
  color0 <- "steelblue"
  color1 <- "deepskyblue3"
  # Plot
  ggplot(eigen_spectrum) + 
    #geom_circle(mapping = aes(x0 = x0, y0 = y0, r = r), data = circle, color = color0) +
    geom_point(mapping = aes(x = Re, y = Im, color = Re), alpha = 0.75) +
    scale_color_continuous(type = "viridis") +
    theme(legend.position = "none") +
    labs(x = "Re", y = "Im", title = paste("Spectrum of a ",mat_str,sep = "")) #+
    #xlim(x_range) + 
    #ylim(-r,r) + 
    #coord_fixed(ratio = 1)
}

#=================================================================================#
#                               HELPER FUNCTIONS
#=================================================================================#

# Read in the eigenvalue in the Kth row from a eigenvalue array and return a numerical
.read_eigenvalue <- function(spectrum, K){complex(real = spectrum[K,1], imaginary = spectrum[K,2])}
