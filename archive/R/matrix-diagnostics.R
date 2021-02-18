#==========================================================================#
#                         RANDOM MATRIX DIAGNOSTICS 
#==========================================================================#

# Visualize the entries of a matrix as a histrogram
visualize_entries <- function(P){
  ggplot() + geom_histogram(data = data.frame(x = as.vector(P)), aes(x = x))
}

# returns proportion of positive entries of any matrix P
pos_entries <- function(P){
  pos_entries <- length(matrix(P[P[,] > 0], nrow = 1))
  pos_entries/(length(P))   
}

#========================================================================#
#                    NORMAL RANDOM MATRIX DIAGNOSTICS
#========================================================================#

normal_params <- function(entries){
  print(paste("Mean: ",round(mean(entries),3),sep=""))
  print(paste("Standard Deviation: ",round(sqrt(var(entries)),3),sep=""))
}

# For a given normal matrix, visualize its entries as a histogram
visualize_normal_entries <- function(P){
  # Vectorize the matrix into a row vector of its entries
  elements_P <- data.frame(x = as.vector(P))
  reals_P <- Re(elements_P$x)
  imags_P <- Im(elements_P$x)
  entries_P <- data.frame(x = c(reals_P, imags_P))
  mu <- mean(entries_P$x)
  print(mu)
  sd <- sqrt(var(entries_P$x))
  print(sd)
  # Get theoretical distribution function
  normal_density <- function(x){dnorm(x, mean = mu, sd = sd)}
  # Plot the histogram of its entries
  entries_hist <- ggplot(data = entries_P, mapping = aes(x)) + 
    geom_histogram(bins = 20, aes(y = stat(density))) +
    stat_function(fun = normal_density)
  # Return plot
  entries_hist
}

#=========================================================================#
#                 STOCHASTIC RANDOM MATRIX DIAGNOSTICS 
#=========================================================================#

# Check if a matrix is stochastic
is_row_stochastic <- function(P){
  row_is_stoch <- rep(F, nrow(P))
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    row_is_stoch[i] <- (row_sum == 1)
  }
  !(F %in% row_is_stoch)
}

# Print out the sums of each of the rows (to check row-stochasticity)
row_sums <- function(P){
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    print(paste("Row ",i,": ",(row_sum),sep=""))
  }
}

# Print out the sums of each of the rows (to check row-stochasticity)
col_sums <- function(P){
  for(i in 1:ncol(P)){
    col_sum <- sum(P[,i])
    print(paste("Col ",i,": ",(col_sum),sep=""))
  }
}