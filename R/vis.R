
#=================================================================================#
#                       RATIO VISUALIZATION FUNCTIONS
#=================================================================================#

# Gives a scatterplot of the ratio entries of the consecutive ratio sequence over time

ratio_scatter <- function(data, n1, mean = 0, range = c(-10,10)){
  ggplot() + 
    geom_scatter(data = data, mapping = aes_string(x= paste("r_x",n1,sep="")), fill = "deepskyblue3") + xlim(range) + 
    geom_vline(xintercept = mean, color = "blue") + scale_fill_discrete(c("")) +
    labs(title = "Distribution of Ratios from the Consecutive Ratio Sequence")
}

# Gives a distribution of the ratio entries of the consecutive ratio sequence
ratio_hist <- function(data, n1, mean = 0, range = c(-10,10)){
  ggplot() + 
    geom_histogram(data = data, mapping = aes_string(x= paste("r_x",n1,sep="")), fill = "deepskyblue3") + xlim(range) + 
    geom_vline(xintercept = mean, color = "blue") + scale_fill_discrete(c("")) +
    labs(title = "Distribution of Ratios from the Consecutive Ratio Sequence")
}

#=================================================================================#
#                               PLOTTING FUNCTIONS
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