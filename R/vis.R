
#=================================================================================#
#                         MIXTIME VISUALIZATION FUNCTIONS
#=================================================================================#

mixtime_histogram <- function(batch, bins = NA, mat = ""){
  mixtimes <- batch$mixtime
  mixtimes <- mixtimes[!is.na(mixtimes)]
  if(is.na(bins)){bins <- range(mixtimes)[2] - range(mixtimes)[1] + 1}
  ggplot(data = batch) + 
    geom_histogram(mapping = aes(x = mixtime), fill = "deepskyblue3", bins = bins) + 
    labs(title = paste("Mixing Time Distribution for a ",mat,"Matrix",sep=""), y = "")
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

scatter_6d <- function(batch){
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

# Gives a scatterplot of the ratio entries of the consecutive ratio sequence over time

ratios_scatterplot <- function(data, n1, mean = 0, range = c(-10,10)){
  ggplot() + 
    geom_scatter(data = data, mapping = aes_string(x= paste("r_x",n1,sep="")), fill = "deepskyblue3") + xlim(range) + 
    geom_vline(xintercept = mean, color = "blue") + scale_fill_discrete(c("")) +
    labs(title = "Distribution of Ratios from the Consecutive Ratio Sequence")
}

# Gives a distribution of the ratio entries of the consecutive ratio sequence
ratios_histogram <- function(ratios, bins = NA){
  mean <- mean(ratios$ratio)
  sd <- var(ratios$ratio)
  range <- c(mean - 3*sd, mean + 3*sd)
  if(is.na(bins)){bins <- 30*range[2]}
  ggplot() + 
    geom_histogram(data = ratios, mapping = aes(x=ratio), fill = "deepskyblue3", bins = bins) +
    geom_vline(xintercept = mean, color = "blue") + scale_fill_discrete(c("")) +
    labs(title = "Distribution of Ratios from the Consecutive Ratio Sequence") +
    xlim(range)
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
  if(nrow(P) == ncol(P)){
    array <- spectrum(P)
  } else{
    array <- P 
  }
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

