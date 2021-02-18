
# This script includes functions that help extract and plot eigenvalues of matrices.


#=================================================================================#
#                              SPECTRUM FUNCTIONS
#=================================================================================#

# Returns a tidied dataframe of the eigenvalues of a matrix ensemble
ensemble_spectrum <- function(ensemble, indexed = FALSE){
  K <- length(ensemble) # Get size of ensemble
  spectra <- spectrum(ensemble[[1]], indexed) # Initialize the spectra stack by evaluating the initial matrix
  # Evaluate the rest of the spectra for the matrices in the ensemble
  for(i in 2:K){
    curr <- spectrum(ensemble[[i]], indexed)
    spectra <- rbind(curr,spectra)
  }
  spectra # Return the spectra for this ensemble
}

# Returns a tidied dataframe of the eigenvalues of a random matrix
spectrum <- function(P, indexed = TRUE){
  M <- nrow(P) # Obtain dimension
  eigen_array <- data.frame(eigen(P)$values) # Get eigenvalues
  eigenvalues <- matrix(rep(NA, 2*M), ncol = 2) # Create matrix to hold eigenvalues
  colnames(eigenvalues) <- c("Re", "Im") # Rename columns
  # Add the components to the array
  for(i in 1:M){
    curr <- eigen_array[i,]
    eigenvalues[i, ] <- c(round(Re(curr), 5), round(Im(curr), 5)) # Round to nearest 5 decimalds
  }
  # Index the eigenvalues
  if(indexed){eigenvalues <- cbind(eigenvalues, data.frame(eigen_index = 1:nrow(eigenvalues)))}
  data.frame(eigenvalues)
}

#=================================================================================#
#                         SPECTRUM VISUALIZATION FUNCTIONS
#=================================================================================#

# Plots the eigenvalues of a given matrix P
spectrum_plot <- function(P, mat_str = ""){
  # See if we have a ensemble of matrices or a single matrix
  not_ensemble <- (nrow(P) == ncol(P))
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
