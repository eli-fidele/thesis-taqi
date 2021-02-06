
#=================================================================================#
#                                 BATCH GENERATION
#=================================================================================#

# Generate and evolve a batch of points for a given matrix P 
run_batch <- function(P, B = 50, lambda = 1, steps = 25, with_ratios = TRUE, final_time = FALSE){
  M <- ncol(P)
  # Make the batch
  batch <- make_batch(M, B, lambda)
  # Evolve the batch and return it
  evolve_batch(P, batch, steps, with_ratios = with_ratios, final_time = final_time)
}

#=================================================================================#
#                       ELEMENTARY BATCH SIMULATION FUNCTIONS
#=================================================================================#

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

# Evolve each element of the batch by a given number of steps and return the evolved stack of arrays
evolve_batch <- function(P, batch, steps, burn_in = 1, with_ratios = TRUE, final_time = FALSE){
  evolved_stack <- evolve(batch[1,], P, steps, burn_in) # Initialize by append first batch element's evolution array
  B <- nrow(batch) # Get number of batch elements
  for(i in 2:B){ 
    #evol <-  # Obtain evolution array of current element of the batch 
    evolved_stack <- rbind(evolved_stack, evolve(batch[i,], P, steps, burn_in)) # Recursively row bind the stack
  }
  rownames(evolved_stack) <- 1:nrow(evolved_stack) # Standardize row names
  evolved_stack <- indexed_batch(evolved_stack, steps) # Index the batch elements 
  if(with_ratios){evolved_stack <- append_ratios(evolved_stack)} # Append ratios if prompted
  if(final_time){evolved_stack <- by_time(evolved_stack, at_time = steps)} # If prompted, return the array at just the final time 
  # Return the stack
  evolved_stack
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
    evolved_row <- as.numeric(v) %*% matrix.power(P,burn_in*i)
    # Evolve the vector, then append the time index to it
    vals[i, ] <- cbind(evolved_row, i)
  }
  # Add the initial batch element to the beginning of the array
  vals <- rbind(c(as.numeric(v),0),vals)
  # Rename steps column to 'time'
  colnames(vals)[ncol(vals)] <- "time"
  vals
}

#=================================================================================#
#                            INDEXING HELPER FUNCTIONS
#=================================================================================#

# This method adds an index to clarify which rows belong to which batch element's evolution array it belongs
indexed_batch <- function(evolved_batch, steps){
  # See if the batch is indexed already
  if(!("element_index" %in% colnames(evolved_batch))){
  # Create the element index column
  element_index <- rep(NA, nrow(evolved_batch))
  # Index the elements in the entire batch using the floor function and return the binded dataframe
  for(i in 1:nrow(evolved_batch)){element_index[i] <- 1 + floor((i-1)/(steps+1))}
  return(data.frame(cbind(evolved_batch,element_index)))
  } # If the batch was indexed already, just return it
  else{return(evolved_batch)}
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
