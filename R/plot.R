
#=================================================================================#
#                         DISPERSION VISUALIZATION FUNCTIONS
#=================================================================================#

#' @title Visualize a plot of the eigenvalue difference spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a histogram of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#' @param metric a string denoting the eigenvalue dispersion metric (column of a dispersion object) to use
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#' @param bins (optional) a string argument of the class of the matrix to label the plot title.
#'
#' Remove "no global binding for variable found" errors from using dplyr verbs by adding .data$ pronoun.
#' @importFrom rlang .data
#'
#' @return A ggplot object containing a histogram of the matrix/matrix ensemble's eigenvalue spacings.
#' @examples
#' # Eigenvalue spacings plot of a normal matrix
#' P <- RM_norm(N = 5)
#' #dispersion.histogram(P)
#'
#' # Eigenvalue spacings plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #dispersion.histogram(Q, mat_str = "Beta")
#'
#' # Eigenvalue spacings plot of an ensemble of normal matrices
#' # ensemble <- RME_norm(N = 3, size = 10)
#' # dispersion.histogram(ensemble)
#'
dispersion.histogram <- function(array, metric = NA, ..., bins = 100){
  valid_schemes <- c("id_diff","id_diff_norm","abs_diff") # Valid schemes for printing if user is unaware of options
  if(class(metric) == "logical"){stop("Please input a valid dispersion metric. Try one of the following: ",paste(valid_schemes, collapse = ", "),".", sep = "")}
  # Process spectrum of the matrix/ensemble
  if(class(array) == "list" || class(array) == "matrix"){ disps_df <- dispersion(array, ...) }
  else{disps_df <- array} # Otherwise, the array is a precomputed dispersion dataframe
  num_entries <- nrow(disps_df) # Get number of entries
  # Plot parameters
  color0 <- "darkorchid4"
  # Return plot
  ggplot(data = disps_df, mapping = aes_string(x = metric)) +
    geom_histogram(mapping = aes(y = stat(count / num_entries)), bins = bins) +
    labs(title = "Distribution of Eigenvalue Spacings", y = "Probability")
}

#' @title Visualize a plot of the eigenvalue difference spectrum of a matrix or ensemble of matrices.
#'
#' @description Returns a scatterplot of the eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#' @param metric a string denoting the eigenvalue dispersion metric (column of a dispersion object) to use
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#'
#' Remove "no global binding for variable found" errors from using dplyr verbs by adding .data$ pronoun.
#' @importFrom rlang .data
#'
#' @return A ggplot object containing a scatterplot of the matrix/matrix ensemble's eigenvalue spacings.
#' @examples
#' # Eigenvalue spacings plot of a normal matrix
#' P <- RM_norm(N = 5)
#' #dispersion.scatterplot(P)
#'
#' # Eigenvalue spacings plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #dispersion.scatterplot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spacings plot of an ensemble of normal matrices
#' # ensemble <- RME_norm(N = 3, size = 10)
#' # dispersion.scatterplot(ensemble)
#'
dispersion.scatterplot <- function(array, metric = "id_diff_norm", pairs = NA, ...){
  # Process dispersion of the matrix/ensemble; if array is a dispersion data frame, copy it.
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, pairs, ...)}
  else{disps_df <- array} # Otherwise, the array is a precomputed dispersion dataframe
  # Parse plotting aesthetics from pairs. diff_ij is more useful unless pairs = "consecutive" or "largest", where j is better.
  if(pairs %in% c("consecutive","largest")){order_stat <- "j"} else{order_stat <- "diff_ij"}
  
  # Plot parameters
  color0 <- "darkorchid4"
  real_valued <- T
  # Scatterplot of dispersion metric
  if(real_valued){  
    disps_df %>%
      ggplot(mapping = aes_string(x = metric, y = order_stat, color = order_stat)) +
      geom_point() +
      scale_color_continuous(type = "viridis") +
      labs(title = "Distribution of Eigenvalue Spacings by Order Statistic")
  } else{
    resolved <- .resolve_eigenvalues(disps_df)
    resolved %>%
      ggplot() +
      geom_point()
  }
}


#' @title Visualize a plot of the variances of the eigenvalue dispersions by order difference level given a matrix or an ensemble.
#'
#' @description Returns a variance scatterplot of classes/levels of eigenvalue spacings of a random matrix or ensemble.
#'
#' @inheritParams dispersion
#' @param metric a string denoting the eigenvalue dispersion metric (column of a dispersion object) to use
#' @param ... any default-valued parameters taken as arguments by spectrum(array, ...)
#'
#' Remove "no global binding for variable found" errors from using dplyr verbs by adding .data$ pronoun.
#' @importFrom rlang .data
#'
#' @return A ggplot object containing a scatterplot of the matrix/matrix ensemble's eigenvalue spacings.
#' @examples
#' # Eigenvalue spacings plot of a normal matrix
#' P <- RM_norm(N = 5)
#' #.dispersion.varplot(P)
#'
#' # Eigenvalue spacings plot of a beta matrix
#' Q <- RM_beta(N = 4, beta = 2)
#' #.dispersion.varplot(Q, mat_str = "Beta")
#'
#' # Eigenvalue spacings plot of an ensemble of normal matrices
#' # ensemble <- RME_norm(N = 3, size = 10)
#' # .dispersion.varplot(ensemble)
#'
.dispersion.varplot <- function(array, metric, ...){
  # Process dispersion of the matrix/ensemble; if array is a dispersion data frame, copy it.
  if(class(array) == "list" || class(array) == "matrix"){disps_df <- dispersion(array, ...)}
  else{disps_df <- array}
  # Plot parameters
  color0 <- "darkorchid4"
  # Get variances by level
  disps_df %>%
    group_by(.data$diff_ij) %>% # Use j for consecutive or largest?
    summarize(Var_Disp = var(.data$Dispersion), size = n()) %>%
    ggplot(mapping = aes(x = .data$diff_ij, y = .data$Var_Disp, color = .data$Var_Disp, size = .data$size)) +
    geom_point() +
    scale_color_continuous(type = "viridis") +
    #scale_size_manual(values = c(1,2)) +
    labs(title = "Variance of Eigenvalue Spacings by Ranking Difference Class", x = "Ranking Difference", y = "Variance")
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
  if(class(array) == "list" || class(array) == "matrix"){
    array_spectrum <- spectrum(array, ...)
  }
  # Else, the array is a precomputed spectrum (avoid computational waste for multiple visualizations)
  else{
    array_spectrum <- array
  }
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
  num_entries <- nrow(array) # Get number of entries to normalize
  if(class(component) == "logical"){component <- c("Re", "Im")} # Set default to both components
  # Plot lambda function
  component_plot <- function(component){
    # Plot parameters
    component_str <- paste(" (",component,")", sep = "")
    # Plot
    array_spectrum %>%
      ggplot() +
      geom_histogram(mapping = aes_string(x = component), fill = color0, bins = bins) +
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

#=================================================================================#
#                             SPECTRUM ORDER PLOTS
#=================================================================================#

order.scatterplot <- function(spectrum, component){
  spectrum %>%
    ggplot(aes(x = Order, y = {{ component }}, color = Order)) +
    geom_point() +
    scale_color_viridis_c() +
    theme(legend.position = "bottom")
}
order.density <- function(spectrum, component){
  spectrum %>%
    ggplot(mapping = aes(group = Order, x = {{ component }}, color = Order)) + 
    geom_density() +
    scale_color_viridis_c() +
    theme(legend.position = "bottom")
}
order.summary <- function(spectrum, component){
  spectrum %>%
    group_by(Order) %>%
    summarize(
      Mean_Re = mean(Re), Mean_Im = mean(Im), Mean_Norm = mean(Norm),
      Variance_Re = var(Re), Variance_Im = var(Im), Variance_Norm = var(Norm)) %>%
    ggplot(mapping = aes(y = {{ component }}, x = Order, color = Order)) + 
    geom_point() +
    geom_line() +
    scale_color_viridis_c() +
    theme(legend.position = "bottom")
}