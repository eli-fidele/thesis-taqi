
#=================================================================================#
#                                 BATCH GENERATION
#=================================================================================#

# Generate and evolve a batch of points for a given matrix P 
run_batch <- function(P, B = 50, lambda = 1, steps = 25, with_ratios = TRUE){
  M <- ncol(P)
  # Make the batch
  batch <- make_batch(M, B)
  # Evolve the batch
  evolved_batch <- evolve_batch(batch, steps, lambda)
  # Add indexing if prompted
  if(with_ratios){evolved_batch <- append_ratios(evolved_batch)} # Append the intermediate ratios
  evolved_batch
}

#=================================================================================#
#                       ELEMENTARY BATCH SIMULATION FUNCTIONS
#=================================================================================#

# Evolve each element of the batch by a given number of steps and return the evolved stack of arrays
evolve_batch <- function(batch, steps, burn_in = 1){
  evolved_stack <- evolve(batch[1,], P, steps, burn_in) # Initialize by append first batch element's evolution array
  B <- nrow(batch) # Get number of batch elements
  for(i in 2:B){ 
    evol <- evolve(batch[i,], P, steps, burn_in) # Obtain evolution array of current element of the batch 
    evolved_stack <- rbind(evolved_stack, evol) # Recursively row bind the stack
  }
  rownames(evolved_stack) <- 1:nrow(evolved_stack) # Standardize row names
  evolved_stack <- indexed_batch(evolved_stack, steps) # Index the batch elements 
  data.frame(evolved_stack)
}

# Generate a Monte Carlo batch
make_batch <- function(M, B, lambda = 1, complex = FALSE){
  batch <- matrix(rep(NA, B * M), nrow = B)  # create [B x M] batch matrix
  # If prompted, generate complex-valued random elements 
  if(complex){ 
    for(i in 1:B){batch[i,] <- complex(real = runif(n = M, min = -lambda, max = lambda), imaginary = runif(n = M, min = -lambda, max = lambda))}
  } else {
      for(i in 1:B){batch[i,] <- runif(n = M, min = -lambda, max = lambda)} # Otherwise, generate real-valued random elements
    }
  batch <- standardize_colnames(batch) # standardize the column names
  data.frame(batch) # return batch
}

# Evolves an element of a batch by the matrix P and returns the evolved array
evolve <- function(v, P, steps, burn_in = 1){
  M <- ncol(P)
  # Simulate the evolution matrix of a given batch element
  vals <- matrix(rep(NA, M * steps), ncol = M)
  vals <- standardize_colnames(vals) # Standardize column names
  # Add a column to track the steps/time
  vals <- cbind(vals, rep(0, steps+1))
  # Evolve the batch element
  for(i in 1:steps){
    # Evolved vector (v * P^steps)
    evolved_row <- as.numeric(v) %*% matrix.power(P,burn_in*i)
    # Add time index to evolved row
    vals[i, ] <- cbind(evolved_row, i)
  }
  #store the values in a dataframe
  vals <- rbind(c(as.numeric(v),0),vals)
  # Rename steps column to 'time'
  colnames(vals)[ncol(vals)] <- "time"
  vals
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

#=================================================================================#
#                            INDEXING HELPER FUNCTIONS
#=================================================================================#

# This method adds an index to clarify which rows belong to which batch element's evolution array it belongs
indexed_batch <- function(evolved_batch, steps){
  # create index column
  element_index <- rep(NA, nrow(evolved_batch))
  # append indexing by running loops
  for(i in 1:nrow(evolved_batch)){element_index[i] <- 1 + floor((i-1)/(steps+1))}
  # only append index column if it is not there already
  if(!("element_index" %in% colnames(evolved_batch))){evolved_batch <- cbind(evolved_batch,element_index)}
  evolved_batch
}

standardize_colnames <- function(array, prefix = ""){
  # get string vector
  str_vec <- rep(NA, ncol(array))
  # rename the columns
  for(i in 1:ncol(array)){
    str_vec[i] <- paste(prefix,"x",i,sep="")
  }
  colnames(array) <- str_vec
  # return renamed col array
  array
}
