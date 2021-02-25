
# This script includes functions that help extract and plot eigenvalues of matrices.

#=================================================================================#
#                              EIGENVALUE DISPERSION
#=================================================================================#

#' @title Obtain the eigenvalue spacings of a matrix or ensemble of matrices.
#'
#' @description Returns a vector of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalue spacings are to be returned
#' @param components returns the array with resolved real and imaginary components; otherwise returns complex-valued vectors/distances
#' @param norm use the norm metric for eigenvalue spacing; otherwise returns absolute difference metric
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
dispersion <- function(array, components = T, norm = T){
  is_ensemble <- (class(array) == "list") # Infer type of array (matrix or ensemble) then parse accordingly
  # Once type of array is inferred, obtain the dispersion array
  if(!is_ensemble){diffs <- .eigen_deltas(array, norm)} # Array is a matrix
  else{diffs <- purrr::map_dfr(.x = array, .f = .eigen_deltas)} # Array is an ensemble; recursively row binding each matrix's dispersions
  diffs # Return differences
}

# Find the eigenvalue dispersion of a given matrix
.eigen_deltas <- function(P, norm = T){
  eigenvalues <- spectrum(P)
  N <- nrow(P)
  diffs <- rep(0, N*(N-1)/2) # N choose 2 number of possible eigenvalue pair differences
  k <- 0
  # Run over pair combination, asserting i > j to avoid repeats
  for(i in 1:N){
    for(j in 1:N){
      if(i < j){
        k <- k + 1
        if(!norm){diffs[k] <- eigenvalues[i,] - eigenvalues[j,]}
        else{diffs[k] <- abs(eigenvalues[i,] - eigenvalues[j,])}
      }
    }
  }
  data.frame(diff = diffs) # Return eigenvalue differences
}

#=================================================================================#
#                         DISPERSION VISUALIZATION FUNCTIONS
#=================================================================================#

#' @title Visualize a plot of the eigenvalue difference spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a scatterplot of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @param array a square matrix or matrix ensemble whose eigenvalues spacings are to be plotted
#' @param bins (optional) a string argument of the class of the matrix to label the plot title.
#'
#' @return A ggplot object containing a scatterplot of the matrix/matrix ensemble's eigenvalue spacings.
#' @examples
#' # Eigenvalue spacings plot of a normal matrix
#' P <- RM_norm(N = 5)
#' #dispersion.plot(P)
#'
#' # Eigenvalue spacings plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #dispersion.plot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spacings plot of an ensemble of normal matrices
#' # ensemble <- RME_norm(N = 3, size = 10)
#' # dispersion.plot(ensemble)
#'
dispersion.plot <- function(array, bins = 100){
  entries <- dispersion(array)
  num_entries <- length(entries) # Get number of entries
  # Plot parameters
  color0 <- "darkorchid4"
  # Return plot
  ggplot(data = entries, aes(x = diff)) +
    geom_histogram(mapping = aes(y = stat(count / num_entries)), fill = color0, bins = bins)+
    scale_fill_discrete(c("")) +
    labs(title = "Distribution of Eigenvalue Spacings", y = "Probability")
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
spectrum <- function(array, components = T, largest = F, smallest = F){
  # Infer type of array (matrix or ensemble) then parse accordingly.
  is_ensemble <- (class(array) == "list")
  round <- 5
  # One type of array is inferred, obtain the eigenvalue array
  if(!is_ensemble){
    P <- array
    eigen_array <- data.frame(eigen(P)$values) # Get eigenvalues of matrix
    spectrum_row <- function(i, tbl){c(round(Re(tbl[i,]), round), round(Im(tbl[i,]), round), abs(tbl[i,]), i)}
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
  if(smallest){eigenvalues <- .specSMALLEST(spectrum)}
  if(largest){eigenvalues <- .specLARGEST(spectrum)}
  #if(class(order) != NULL){eigenvalues <- eigenvalues %>% filter(Order %in% c(order))}
  # Return spectrum of eigenvalues
  eigenvalues
}

# Returns smallest eigenvalues of the matrix
.specSMALLEST <- function(spectrum){spectrum[which(spectrum$Order == which.max(spectrum$Order)),]}

# Returns largest eigenvalues of the matrix
.specLARGEST <- function(spectrum){spectrum[which(spectrum$Order == 1),]}

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
#' # Eigenvalue spectrum plot of a matrix
#' P <- RM_norm(N = 5)
#' spectrum.plot(P)
#'
#' # Labelled spectrum plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' spectrum.plot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spectra plot of an ensemble of normal matrices
#' ensemble <- RME_norm(N = 3, size = 10)
#' spectrum.plot(ensemble)
#'
spectrum.plot <- function(array, mat_str = ""){
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){eigen_spectrum <- spectrum(array)}
  else{eigen_spectrum <- array}
  if(mat_str != ""){pre_space <- " "} else{pre_space <- ""} # Format without given name
  # Infer plot title string from which type of array (matrix/ensemble)
  is_mat <- class(array) == "matrix"
  if(is_mat){plot_str <- paste(pre_space, mat_str," Matrix", sep = "", collapse = "")}
  else{plot_str <- paste(pre_space, mat_str," Matrix Ensemble", sep = "", collapse = "")}
  # Plot parameters
  #r <- 1
  #x_window <- 0.5
  #x_range <- c(-(r + x_window), (r + x_window)) # Widen the width of the plot
  #circle <- data.frame(x0 = 0, y0 = 0, r = r)
  # Color plot parameters
  #color0 <- "steelblue"
  #color1 <- "deepskyblue3"
  #panel0 <- "lightblue"
  #panel1 <- "lightskyblue1"
  order <- eigen_spectrum[["Order"]]
  # Plot
  ggplot2::ggplot(eigen_spectrum) +
    #geom_circle(mapping = aes(x0 = x0, y0 = y0, r = r), data = circle, color = color0) +
    geom_point(mapping = aes(x = Re, y = Im, color = order), alpha = 0.75) +
    scale_color_continuous(type = "viridis") +
    labs(x = "Re", y = "Im", title = paste("Spectrum of a",plot_str,sep = "")) #+
    #theme(legend.position = "none") +
    #theme(
    #  panel.background = element_rect(fill = panel0,
    #                                  colour = panel0,
    #                                  size = 0.5, linetype = "solid"),
    #  panel.grid.major = element_line(size = 0.5, linetype = 'solid',
    #                                  colour = panel1),
    #  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
    #                                  colour = panel1))#+
    #xlim(x_range) +
    #ylim(-r,r) +
    #coord_fixed(ratio = 1)
}

