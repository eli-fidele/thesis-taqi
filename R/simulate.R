
#=================================================================================#
#                                 BATCH GENERATION
#=================================================================================#

# Simulation to determine the eigenvalue mixing time of a random matrix.
mixtime_sim <- function(P, batch_size, steps, epsilon = 0.1){
  B <- batch_size
  init_sim <- initial_sim(P, B, steps) # Initialize the simulation
  # Extract the arrays
  batch <- init_sim[[1]]
  evolved_batch <- init_sim[[2]]
  # Perform the mixing time analysis
  evolved_batch <- eigen_classify(evolved_batch, P, epsilon) # Classify candidate eigenvectors
  eigen_index <- by_time(evolved_batch, at_time = steps)$eigen_index # Get the array at the final time
  batch <- cbind(batch, eigen_index) # Add the eigen_index of the batch elements to the array
  batch <- eigen_mixtime(evolved_batch, batch) # Add the mixing times to the batch element array
  # Return list of arrays
  list(batch, evolved_batch)
}

# Generate and evolve a batch of points for a given random matrix P. 
# This function is a basic "initial" simulation. Other simulation functions will utilize this function 
initial_sim <- function(P, B, steps){
  # Make the batch
  batch <- make_batch(M = ncol(P), B)
  # Evolve the batch and return it
  evolved_batch <- evolve_batch(batch, P, steps, ratios = TRUE)
  list(batch, evolved_batch)
}

#=================================================================================#
#                       ELEMENTARY BATCH SIMULATION FUNCTIONS
#=================================================================================#

# Generate a Monte Carlo batch
make_batch <- function(M, B, lambda = 1, complex = FALSE){
  batch <- matrix(rep(NA, B * M), nrow = B)  # create [B x M] batch matrix
  if(complex){ # If prompted, generate complex-valued random elements 
    for(i in 1:B){batch[i,] <- complex(real = runif(n = M, min = -lambda, max = lambda), imaginary = runif(n = M, min = -lambda, max = lambda))}
  } else {
    for(i in 1:B){batch[i,] <- runif(n = M, min = -lambda, max = lambda)} # Otherwise, generate real-valued random elements
  }
  batch <- standardize_colnames(batch) # standardize the column names
  data.frame(batch) # return batch
}

# Evolve each element of the batch by a given number of steps and return the evolved stack of arrays
evolve_batch <- function(batch, P, steps, burn_in = 1, ratios = TRUE){
  B <- nrow(batch) # Get number of batch elements
  evolved_stack <- evolve(batch[1,], P, steps, burn_in) # Initialize by append first batch element's evolution array
  for(i in 2:B){ 
    evol <-  evolve(batch[i,], P, steps, burn_in) # Obtain evolution array of current element of the batch 
    evolved_stack <- rbind(evolved_stack, evol) # Recursively row bind the stack
  }
  rownames(evolved_stack) <- 1:nrow(evolved_stack) # Standardize row names
  evolved_stack <- add_indices(evolved_stack, steps) # Index the batch elements 
  if(ratios){evolved_stack <- append_ratios(evolved_stack)} # Append ratios if prompted
  # Return the stack
  evolved_stack
}

# Evolves an element of a batch by the matrix P and returns the array of the evolution sequence
evolve <- function(v, P, steps, burn_in = 1){
  M <- ncol(P)
  # Simulate the evolution matrix of a given batch element
  seq <- matrix(rep(NA, M * steps), ncol = M)
  seq <- standardize_colnames(seq) # Standardize column names
  seq <- cbind(seq, rep(0, steps+1)) # Add a column to track the steps/time
  # Evolve the batch element
  for(i in 1:steps){
    evolved_row <- as.numeric(v) %*% matrix.power(P,burn_in*i)
    # Evolve the vector, then append the time index to it
    seq[i, ] <- cbind(evolved_row, i)
  }
  # Add the initial batch element to the beginning of the array
  seq <- rbind(c(as.numeric(v),0),seq)
  # Rename steps column to 'time'
  colnames(seq)[ncol(seq)] <- "time"
  seq
}

#=================================================================================#
#                            NAMING/INDEXING HELPER FUNCTIONS
#=================================================================================#

# This method adds an index to clarify which rows belong to which batch element's evolution array it belongs
add_indices <- function(evolved_batch, steps){
  # Create the element index column
  element_index <- rep(NA, nrow(evolved_batch))
  # Index the elements in the entire batch using the floor function and return the binded dataframe
  for(i in 1:nrow(evolved_batch)){element_index[i] <- 1 + floor((i-1)/(steps+1))}
  evolved_batch <- cbind(evolved_batch,element_index)
  data.frame(evolved_batch)  # Return the batch
}

# Standardizes an array with n dimensions to have column names of the form "prefix-x(i)" for i = 1,...,n
standardize_colnames <- function(array, prefix = ""){
  # Initialize a string vector for the column names
  new_colnames <- rep("", ncol(array))
  # Generate the standardized column names, with prefix if added
  for(i in 1:ncol(array)){new_colnames[i] <- paste(prefix,"x",i,sep="")}
  colnames(array) <- new_colnames
  array # Return array
}
