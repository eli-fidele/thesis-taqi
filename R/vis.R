
#=================================================================================#
#                         MIXTIME VISUALIZATION FUNCTIONS
#=================================================================================#

# Returns a histogram of the mixing times of a batch array
mixtime_histogram <- function(batch, bins = NA, mat = ""){
  mixtimes <- batch$mixtime # Get vector of mixtimes
  mixtimes <- mixtimes[!is.na(mixtimes)] # Remove NAs
  if(class(bins) == "logical"){bins <- range(mixtimes)[2] - range(mixtimes)[1] + 1}
  # Return plot
  color0 <- "deeporchid2"
  ggplot(data = batch, mapping = aes(x = mixtime)) + 
    geom_histogram(, fill = color0, bins = bins) + 
    labs(title = paste("Mixing Time Distribution for a ",mat,"Random Matrix",sep=""), y = "")
}

mixtime4d <- function(batch, grid = T){
  plot1 <- ggplot() + geom_point(data = batch, mapping = aes(x = x1, y = mixtime, color = mixtime))
  plot2 <- ggplot() + geom_point(data = batch, mapping = aes(x = x2, y = mixtime, color = mixtime))
  plot3 <- ggplot() + geom_point(data = batch, mapping = aes(x = x3, y = mixtime, color = mixtime))
  plot4 <- ggplot() + geom_point(data = batch, mapping = aes(x = x4, y = mixtime, color = mixtime))
  if(grid){
    grid.arrange(plot1,plot2,plot3,plot4)
  } else{
    list(plot1,plot2,plot3,plot4)
  }
}

mixtime6d <- function(batch){
  plot1 <- ggplot() + geom_point(data = batch, mapping = aes(x = x1, y = x2, color = mixtime))
  plot2 <- ggplot() + geom_point(data = batch, mapping = aes(x = x2, y = x3, color = mixtime))
  plot3 <- ggplot() + geom_point(data = batch, mapping = aes(x = x3, y = x4, color = mixtime))
  plot4 <- ggplot() + geom_point(data = batch, mapping = aes(x = x4, y = x2, color = mixtime))
  plot5 <- ggplot() + geom_point(data = batch, mapping = aes(x = x4, y = x3, color = mixtime))
  plot6 <- ggplot() + geom_point(data = batch, mapping = aes(x = x4, y = x1, color = mixtime))
  list(plot1,plot2,plot3,plot4,plot5,plot6)
}


#=================================================================================#
#                         RATIO VISUALIZATION FUNCTIONS
#=================================================================================#

# Gives a variance scatterplot of the ratio entries over time
variance_scatterplot <- function(evolved_batch, at_time = NA, log = T){
  max_time <- max(evolved_batch$time) # Get the maximum time for the evolved batch
  # Initialize the time range to default
  if(class(at_time) == "logical"){at_time <- 2:max_time}
  # Get variances
  variances <- variance_by_time(evolved_batch, at_time, log = T)
  # Setup and return plot
  color0 <- "darkorchid3"
  if(log){plot_str <- "Log-"} else{plot_str <- ""}
  ggplot(data = data.frame(time = at_time, variance = variances),
         mapping = aes(x = time, y = variance)) + 
    geom_point(color = color0) +
    geom_line(color = color0) +
    labs(title = paste(plot_str,"Variance of the Ratio Entries by Matrix Power",sep=""),
         y = paste(plot_str,"Variance",sep=""))
}

# Gives a distribution of the ratio entries of the consecutive ratio sequence
ratios_histogram <- function(ratios, at_time = NA, log = T, alpha = 0.99, bins = 200){
  max_time <- max(evolved_batch$time)
  # Initialize the time range to default
  if(class(at_time) == "logical"){at_time <- 2:max_time}
  ratios <- data.frame(ratio = ratios_by_time(evolved_batch, at_time, log = T))
  # Return plot
  color0 <- "darkorchid3"
  if(log){plot_str <- "Log-"} else{plot_str <- ""}
  ggplot(data = ratios, mapping = aes(x=ratio)) + 
    geom_histogram(fill = color0, bins = bins) +
    scale_fill_discrete(c("")) +
    labs(title = paste("Distribution of ",plot_str,
                       "Ratios from the Consecutive Ratio Sequence",sep="")) +
    xlim(quantile(ratios$ratio, probs = quantiles_alpha(alpha))) # 
}

# Create a vector of quantiles containing alpha % of the data
quantiles_alpha <- function(alpha){
  remaining <- (1 - alpha)/2
  c(remaining, alpha+remaining)
}

# Gives a scatterplot of the ratio entries of the consecutive ratio sequence over time
ratios_scatterplot <- function(data, n1, mean = 0, range = c(-10,10)){
  ggplot() + 
    geom_scatter(data = data, mapping = aes_string(x= paste("r_x",n1,sep="")), fill = "deepskyblue3") + xlim(range) + 
    geom_vline(xintercept = mean, color = "blue") + scale_fill_discrete(c("")) +
    labs(title = "Distribution of Ratios from the Consecutive Ratio Sequence")
}


#=================================================================================#
#                           BATCH VISUALIZATION FUNCTIONS
#=================================================================================#

# Plots the evolution arrays of a 3D evolved batch
batch_3d_customplot <- function(batch_data,n1,n2,n3,mat_str=""){
  plot_empty <- ggplot() 
  plot_12 <- batch_2d_customplot(batch_data, n1, n2,mat_str)
  plot_23 <- batch_2d_customplot(batch_data, n2, n3,mat_str)
  plot_13 <- batch_2d_customplot(batch_data, n1, n3,mat_str)
  grid.arrange(plot_empty,plot_12,plot_23,plot_13, ncol = 2)
}

# Plots the evolution arrays of a 3D evolved batch
batch_3d_plot <- function(batch_data,mat_str=""){
  plot_empty <- ggplot() 
  plot_12 <- batch_2d_customplot(batch_data, 1, 2,mat_str)
  plot_23 <- batch_2d_customplot(batch_data, 2, 3,mat_str)
  plot_13 <- batch_2d_customplot(batch_data, 1, 3,mat_str)
  grid.arrange(plot_empty,plot_12,plot_23,plot_13, ncol = 2)
}

# Plots the evolution arrays of a 2D evolved batch given two particular dimensions
batch_2d_customplot <- function(batch_data, n1, n2, mat_str = ""){
  ggplot(batch_data, mapping = aes(color = as.factor(element_index))) + 
    geom_point(mapping = aes_string(x = paste("x",n1,sep=""), y = paste("x",n2,sep=""))) +
    theme(legend.position = "none") +
    labs(title = paste("Evolution of a Markov Chain",mat_str))
}

# Plots the evolution arrays of a 2D evolved batch
batch_2d_plot <- function(batch_data, mat_str = ""){
  ggplot(batch_data, mapping = aes(x = x1, y = x2, color = as.factor(element_index))) + 
    geom_point() +
    theme(legend.position = "none") +
    labs(title = paste("Evolution of a Markov Chain",mat_str))
}

#=================================================================================#
#                         SPECTRUM VISUALIZATION FUNCTIONS
#=================================================================================#

# Plots the eigenvalues of a given matrix P
spectrum_plot <- function(P, mat_type=""){
  # Check if we have a stack of matrices or singular matrix
  if(nrow(P) == ncol(P)){array <- spectrum(P)} else{array <- P}
  # Plot parameters
  r <- 1
  ep <- 0.5
  # Plot
  ggplot(array) + 
    geom_point(aes(x = Re, y = Im), color = "deepskyblue3") + 
    labs(x = "Re", y = "Im", title = paste("Spectrum of an ",mat_type,"Ensemble",sep = "")) +
    xlim(-(r+ep),(r+ep)) + ylim(-r,r) +
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r), color = "steelblue") +
    coord_fixed(ratio = 1) +
    theme(legend.position = "none")
}

