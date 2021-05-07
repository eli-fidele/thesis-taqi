

#=================================================================================#
#                       ELEMENTARY BATCH SIMULATION FUNCTIONS
#=================================================================================#

# Generate a Monte Carlo batch of uniform points in an M-hypercube, or a random set of initial probability distributions
generate_batch <- function(N, batch_size, lambda = 1, complex = FALSE, stoch = FALSE){
  B <- batch_size # Rename variable for brevity
  # Create random initial probability distributions if matrix is stochastic
  if(stoch){batch <- do.call("rbind", lapply(X = rep(N, B), FUN = .stoch_row))}
  # Otherwise, generate real-valued uniformly random elements
  else{batch <- do.call("rbind", lapply(X = rep(N, B), FUN = function(N, lambda){runif(N, -lambda, lambda)}, lambda = lambda))} 
  if(complex){batch <- batch + 1i * generate_batch(N, B, lambda, stoch = stoch)} # Add complex component if prompted
  batch <- .standardize_colnames(batch) # Standardize the column names
  data.frame(batch) # Return batch
}

# Evolve each element of the batch by a given number of steps and return the evolved stack of arrays
evolve_batch <- function(batch, P, steps){
  # Get number of batch elements
  B <- nrow(batch)
  .evolve_element <- function(i, batch, P, steps){ data.frame(evolve(v = batch[i,], P, steps)) }
  # Iteratively rowbind the evolved row at powers k = 1,...,steps for all the batch elements
  evolved_stack <- purrr::map_dfr(1:B, .evolve_element, batch, P, steps)
  # Preprocess for return
  evolved_stack <- .add_indices(evolved_stack, steps) # Index the batch elements 
  evolved_stack <- .append_ratios(evolved_stack) # Append ratios
  rownames(evolved_stack) <- 1:nrow(evolved_stack) # Standardize row names
  # Return the stack
  evolved_stack 
}

# Evolves an element of a batch (v) by the matrix P and returns the array of the evolution sequence
evolve <- function(v, P, steps){
  # Call the matrix over 0:steps, 0 being the original element
  seq <- do.call("rbind",lapply(X = 0:steps, FUN = .MATpower, v = v, P = P)) 
  seq <- .standardize_colnames(seq, time = T) # Standardize column names and labelling the time column
  seq # Return the array 
}

# Returns vector multiplied by k^th power of P and index of k as "time"
.MATpower <- function(k, v, P){c((as.numeric(v) %*% matrix.power(P,k)), k)}


#=================================================================================#
#                                 HELPER FUNCTIONS
#=================================================================================#

# Append ratio of row elements by each step
.append_ratios <- function(evolved_batch){
  # Extract number of batch elements
  B <- max(evolved_batch$element_index)
  ratio_stack <- do.call("rbind",lapply(1:B, .ratios_by_element, evolved_batch)) # Get the element ratios for every element
  ratio_stack <- .standardize_colnames(ratio_stack, prefix = "r_") # Standardize the column names
  # Return evolved batch with the ratios
  cbind(evolved_batch, ratio_stack)
}

# Find the element ratios array between the steps for a given element array for all the times t = 1,...,steps
.ratios_by_element <- function(element_index, evolved_batch){
  curr_element <- by.element(evolved_batch, element_index) # Filter for current element
  # Extract parameters
  N <- ncol(curr_element) - 2 # No. of dimensions = No. of cols - time & element_index columns
  steps <- nrow(curr_element) - 1 # Initial element not counted
  # Helper function extracting the entry-wise ratio of elements for an element's entries at t = i by t = i-1.
  .ROWratio <- function(i, evolved_element, N){evolved_element[i+1, 1:N]/evolved_element[i, 1:N]}
  # Get rows for currently indexed element
  evolved_element <- by.element(evolved_batch, element_index) 
  # Rowbind the ratios at time = 1 to time = steps.
  rbind(rep(NA, N), do.call("rbind",lapply(X = 1:steps, FUN = .ROWratio, evolved_element, N))) # Return element ratio stack
}

# Extract the array for a particular element/a range of elements
by.element <- function(array, index){
  array[which(array$element_index %in% c(index)),]
}

# Extract the array for a particular time/a range of times
by.time <- function(array, at_time){
  array[which(array$time %in% c(at_time)),]
}

#=================================================================================#
#                         NAMING/INDEXING HELPER FUNCTIONS
#=================================================================================#

# This function take the row of a evolved array with the number of steps used. It infers which element it is indexing and adds the time.
# Since each batch element is evolved steps number of times, we use the floor function to map {1,...,steps+1} to {0,...,steps}.
.add_indices <- function(evolved_batch, steps){
  num_states <- nrow(evolved_batch) # Get total number of states
  .STATEidx <- function(i, steps){1 + floor((i-1)/(steps+1))} # Lambda helper function; formula for .add_indices
  element_index <- data.frame(element_index = map_dbl(1:num_states, .f = .STATEidx, steps = steps))
  cbind(evolved_batch,element_index) # Return the indexed evolved batch
}

# Standardizes an array with n dimensions to have column names of the form "prefix-x(i)" for i = 1,...,n
.standardize_colnames <- function(array, time = F, prefix = ""){
  N <- ncol(array) - as.numeric(time) # Get dimension of the matrix
  .DIMstr <- function(i, prefix = ""){paste(prefix,"x",i,sep="")} # Standarized element dimension column name ("xi" for DIM i)
  # Map each column to standarized dimension string and add time column if prompted
  if(time){colnames(array) <- c(purrr::map_chr(1:N, .f = .DIMstr, prefix = prefix), "time")}
  else{colnames(array) <- purrr::map_chr(1:N, .f = .DIMstr, prefix = prefix)}
  array # Return array
}

#=================================================================================#
#                         STOCHASTIC HELPER FUNCTIONS
#=================================================================================#

# Generates stochastic rows of length N
.stoch_row <- function(N){
  # Sample a vector of probabilities
  row <- runif(n = N, min = 0, max = 1)
  # Return the normalized row (sums to one)
  row / sum(row)
}

# Check if a matrix is stochastic
.isStochastic <- function(P){
  row_is_stoch <- rep(F, nrow(P))
  for(i in 1:nrow(P)){
    row_sum <- sum(P[i,])
    row_is_stoch[i] <- (row_sum == 1)
  }
  !(F %in% row_is_stoch)
}