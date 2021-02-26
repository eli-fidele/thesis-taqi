
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
  is_ensemble <- (class(array) == "list") # Infer type of array (matrix or ensemble)
  # Array is a matrix; call function returning dispersion for singleton
  if(!is_ensemble){.dispersion_matrix(array, norm)}
  # Array is an ensemble; recursively row binding each matrix's dispersions
  else{
    pairs <- .unique_pairs(N) # Compute pairs to avoid computational waste and pass as argument
    purrr::map_dfr(.x = array, .f = .dispersion_matrix, norm, components, digits, pairs)
  }
}

# Find the eigenvalue dispersions for a given matrix
.dispersion_matrix <- function(P, norm = T, components = T, digits = 3, pairs = NA){
  #eigenvalues <- .spectrum_matrix(P, components = F, digits) # Get the eigenvalues of a matrix
  eigenvalues <- eigen(P)$values
  N <- nrow(P) # Get matrix dimension
  # If uninitialized for the ensemble, enumerate unique pairs of N eigenvalues
  if(class(pairs) == "logical"){idx_pairs <- .unique_pairs(N)}
  else{idx_pairs <- pairs} # Otherwise, read in pre-computed values 
  # User is requesting a norm function rather than a setting of Euclidean norm
  if(class(norm) != "logical"){
    norm_fn <- function(x){(abs(x))^norm}
    purrr::map2_dfr(idx_pairs[,1], idx_pairs[,2], .resolve_dispersion, eigenvalues, norm_fn, components, digits)
  } else{
    purrr::map2_dfr(idx_pairs[,1], idx_pairs[,2], .resolve_dispersion, eigenvalues, norm, components, digits)
    }
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
  cbind(disp, data.frame(OrderDiff = i - j))
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
  ggplot(data = entries, aes(x = Dispersion)) +
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
  # Infer type of array (matrix or ensemble) then parse accordingly.
  is_ensemble <- (class(array) == "list")
  # One type of array is inferred, obtain the eigenvalue array
  if(!is_ensemble){.spectrum_matrix(array)}
  # Otherwise, recursively get ensemble's spectrum by row binding each matrix's spectrum
  else{purrr::map_dfr(1:length(array), FUN = function(i){.spectrum_matrix}, components, largest, smallest, digits)}
}

# Helper function returning tidied eigenvalue array for a matrix
.spectrum_matrix <- function(P, components = T, largest = F, smallest = F, digits = 3){
  eigenvalues <- eigen(P)$values # Get eigenvalues of matrix P
  # Get largest eigenvalue
  if(largest){.resolve_eigenvalue(order = 1, eigenvalues, components)} 
  # Get smallest eigenvalue
  else if(smallest){.resolve_eigenvalue(order = nrow(P), eigenvalues, components)}
  # Get all the eigenvalues
  purrr::map_dfr(1:nrow(P), .resolve_eigenvalue, eigenvalues, components, digits)
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

