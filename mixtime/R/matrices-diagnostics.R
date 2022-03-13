#==========================================================================#
#                         RANDOM MATRIX DIAGNOSTICS 
#==========================================================================#

# Visualize the entries of a matrix as a histrogram
.visualize_entries <- function(P){
  ggplot() + geom_histogram(data = data.frame(x = as.vector(P)), aes(x = x))
}

#========================================================================#
#                    NORMAL RANDOM MATRIX DIAGNOSTICS
#========================================================================#

.normal_params <- function(entries){
  # Vectorize the matrix into a row vector of its entries
  if(class(entries) == "matrix"){
    elements_P <- data.frame(x = as.vector(P))
    entries <- data.frame(x = c(Re(elements_P$x), Im(elements_P$x)))
    }
  print(paste("Mean: ",round(mean(entries$x),3),sep=""))
  print(paste("Standard Deviation: ",round(sqrt(var(entries)),3),sep=""))
}

# For a given normal matrix, visualize its entries as a histogram
.visualize_normal_entries <- function(P){
  # Vectorize the matrix into a row vector of its entries
  elements_P <- data.frame(x = as.vector(P))
  entries_P <- data.frame(x = c(Re(elements_P$x), Im(elements_P$x)))
  mu <- mean(entries_P$x, na.rm = T)
  sd <- sqrt(var(entries_P$x, na.rm = T))
  .normal_params(entries_P)
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

# Print out the sums of each of the rows (to check row-stochasticity)
.row_sums <- function(P){for(i in 1:nrow(P)){print(paste("Row ",i,": ",(sum(P[i,])),sep=""))}}

# Print out the sums of each of the rows (to check row-stochasticity)
.col_sums <- function(P){for(i in 1:ncol(P)){print(paste("Col ",i,": ",(sum(P[,i])),sep=""))}}

