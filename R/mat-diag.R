#==========================================================================#
#                         RANDOM MATRIX DIAGNOSTICS
#==========================================================================#

# Visualize the entries of a matrix as a histrogram
.visualize_entries <- function(P){
  # Vectorize the matrix into a row vector of its entries
  entries <- as.vector(P)
  if(class(entries) == "complex"){
    entries <- data.frame(entries = c(Re(entries), Im(entries)))
  } else{
    entries <- data.frame(entries = entries)
  }
  # Plot the histogram of its entries
  entries_hist <-
    ggplot(data = entries) +
    geom_histogram(bins = 20, aes(x = entries, y = stat(density)))
  # Return plot
  entries_hist
}

#========================================================================#
#                    NORMAL RANDOM MATRIX DIAGNOSTICS
#========================================================================#

.normal_params <- function(array){
  if(class(array) == "data.frame"){entries <- array[[1]]}
  else{entries <- as.vector(array)}
  mu <- mean(entries)
  sd <- sqrt(var(entries))
  # Vectorize the matrix into a row vector of its entries
  #print(paste("Mean: ",round(mu, 3),sep=""))
  #print(paste("Standard Deviation: ",round(sd, 3),sep=""))
  list(mu, sd)
}

# For a given normal matrix, visualize its entries as a histogram
.visualize_normal_entries <- function(P){
  # Vectorize the matrix into a row vector of its entries
  entries <- as.vector(P)
  if(class(entries) == "complex"){
    entries <- data.frame(entries = c(Re(entries), Im(entries)))
  } else{
    entries <- data.frame(entries = entries)
  }
  # Get parameters
  mu <- mean(entries[[1]], na.rm = T)
  sd <- sqrt(var(entries[[1]], na.rm = T))
  .normal_params(entries) # Print parameters
  # Get theoretical distribution function
  normal_density <- function(x){dnorm(x, mean = mu, sd = sd)}
  # Plot the histogram of its entries
  entries_hist <-
    ggplot(data = entries) +
    geom_histogram(bins = 20, aes(x = entries, y = stat(density))) +
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

