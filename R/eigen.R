
# This script includes functions that help extract and plot eigenvalues of matrices.


#=================================================================================#
#                              SPECTRUM FUNCTIONS
#=================================================================================#

# Returns a tidied dataframe of the eigenvalues of a random matrix
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
  eigenvalues <- matrix(rep(NA,2*M), ncol = 2) # Create matrix to hold eigenvalues
  colnames(eigenvalues) <- c("Re","Im") # Rename columns
  # Add the components to the array
  for(i in 1:M){
    curr <- eigen_array[i,]
    eigenvalues[i, ] <- c(round(Re(curr),5),round(Im(curr),5)) 
  }
  # Index the eigenvalues
  if(indexed){eigenvalues <- cbind(eigenvalues, data.frame(Index = 1:nrow(eigenvalues)))}
  data.frame(eigenvalues)
}

#=================================================================================#
#                               HELPER FUNCTIONS
#=================================================================================#

# Read in the eigenvalue in the Kth row from a eigenvalue array and return a numerical
read_eigenvalue <- function(spectrum, K){complex(real = spectrum[K,1], imaginary = spectrum[K,2])}

#=================================================================================#
#                         SPECTRUM VISUALIZATION FUNCTIONS
#=================================================================================#

# Plots the eigenvalues of a given matrix P
spectrum_plot <- function(P, mat_type=""){
  # Plot parameters
  r <- 1
  ep <- 0.5
  # Check if we have a stack of matrices or singular matrix
  if(nrow(P) == ncol(P)){
    array <- spectrum(P)
  } else{array <- P}
  # Plot
  ggplot(array) + 
    geom_point(aes(x = Re, y = Im), color = "deepskyblue3") + 
    labs(x = "Re", y = "Im", title = paste("Spectrum of an ",mat_type,"Ensemble",sep = "")) +
    xlim(-(r+ep),(r+ep)) + ylim(-r,r) +
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r), color = "steelblue") +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none")
}
