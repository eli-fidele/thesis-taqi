

#=================================================================================#
#                       ELEMENTARY BATCH SIMULATION FUNCTIONS
#=================================================================================#

# Generate a Monte Carlo batch of uniform points in an M-hypercube, or a random set of initial probability distributions
generate_batch <- function(N, batch_size, lambda = 1, complex = FALSE, stoch = FALSE){
  B <- batch_size # Rename variable for aesthetic
  batch <- matrix(rep(NA, B * N), nrow = B)  # create [B x N] batch matrix
  for(i in 1:B){
    # Generate real-valued uniformly random elements unless matrix stochastic
    if(!stoch){batch[i,] <- runif(n = N, min = -lambda, max = lambda)} else{
      batch[i,] <- .stoch_row(N) # Otherwise, create random initial probability distributions
    } 
  }
  if(complex){batch <- batch + 1i * generate_batch(N, B, lambda, stoch = stoch)} # Add complex component if prompted
  batch <- .standardize_colnames(batch) # Standardize the column names
  data.frame(batch) # Return batch
}

# Evolve each element of the batch by a given number of steps and return the evolved stack of arrays
evolve_batch <- function(batch, P, steps){
  B <- nrow(batch) # Get number of batch elements
  evolved_stack <- evolve(batch[1,], P, steps) # Initialize by append first batch element's evolution array
  for(i in 2:B){ 
    evol <-  evolve(batch[i,], P, steps) # Obtain evolution array of current element of the batch 
    evolved_stack <- rbind(evolved_stack, evol) # Recursively row bind the stack
  }
  evolved_stack <- .add_indices(evolved_stack, steps) # Index the batch elements 
  evolved_stack <- .append_ratios(evolved_stack) # Append ratios
  rownames(evolved_stack) <- 1:nrow(evolved_stack) # Standardize row names
  # Return the stack
  evolved_stack
}

# Evolves an element of a batch (v) by the matrix P and returns the array of the evolution sequence
evolve <- function(v, P, steps){
  # Call the matrix over 0:steps, 0 being the original element
  seq <- do.call("rbind",lapply(X = 0:steps, FUN = .mat_power, v = v, P = P)) 
  seq <- .standardize_colnames(seq, time = T) # Standardize column names and labelling the time column
  seq # Return the array 
}


#=================================================================================#
#                                 HELPER FUNCTIONS
#=================================================================================#

# Returns vector multiplied by i^th power of P and index of i as "time"
.mat_power <- function(i, v, P){c((as.numeric(v) %*% matrix.power(P,i)), i)}

# Append ratio of row elements by each step
.append_ratios <- function(evolved_batch){
  # Extract number of batch elements
  B <- max(evolved_batch$element_index)
  # Get the ratios in the array for the first element
  ratio_stack <- .ratios_by_element(evolved_batch, element_index = 1)
  # Repeat for the rest, concatenating by row
  for(i in 2:B){
    curr_ratios <- .ratios_by_element(evolved_batch, element_index = i)
    ratio_stack <- rbind(ratio_stack, curr_ratios)
  }
  # Standardize the column names
  ratio_stack <- .standardize_colnames(ratio_stack, prefix = "r_")
  # Return evolved batch with the ratios
  cbind(evolved_batch, ratio_stack)
}

# Find the ratios between the steps for a given element array
.ratios_by_element <- function(evolved_batch, element_index){
  # Get the array for the current indexed element
  curr_element <- by.element(evolved_batch, element_index)
  # Get number of dimensions and steps
  M <- ncol(curr_element) - 2
  steps <- nrow(curr_element)
  # Initalize the stack
  ratio_stack <- rep(NA, M)
  for(i in 2:steps){
    # Get ratio of rows from current step 
    curr_ratios <- curr_element[i, 1:M]/curr_element[i-1, 1:M]
    # Stack
    ratio_stack <- rbind(ratio_stack, curr_ratios)
  }
  ratio_stack
}

# Extract the array for a particular element/a range of elements
by.element <- function(evolved_batch, index){
  if(class(index) == "numeric"){
    evolved_batch %>% filter(element_index == index)
  } else {
    evolved_batch %>% filter(element_index %in% index)
  }
}

# Extract the array for a particular time/a range of times
by.time <- function(evolved_batch, at_time){
  if(class(at_time) == "numeric"){
    evolved_batch %>% filter(time == at_time)
  } else {
    evolved_batch %>% filter(time %in% at_time)
  }
}


#=================================================================================#
#                         NAMING/INDEXING HELPER FUNCTIONS
#=================================================================================#

# This method adds an index to clarify which rows belong to which batch element's evolution array it belongs
.add_indices <- function(evolved_batch, steps){
  # Create the element index column
  element_index <- rep(NA, nrow(evolved_batch))
  # Index the elements in the entire batch using the floor function and return the binded dataframe
  for(i in 1:nrow(evolved_batch)){element_index[i] <- 1 + floor((i-1)/(steps+1))}
  evolved_batch <- cbind(evolved_batch,element_index)
  data.frame(evolved_batch)  # Return the batch
}

# Standardizes an array with n dimensions to have column names of the form "prefix-x(i)" for i = 1,...,n
.standardize_colnames <- function(array, time = F, prefix = ""){
  N <- ncol(array) - as.numeric(time) # Get dimension of the matrix
  # Initialize a string vector for the column names, adding one colmun if time column is included
  new_colnames <- rep("", N + as.numeric(time))
  # Generate the standardized element column names, with prefix if added
  for(i in 1:N){new_colnames[i] <- paste(prefix,"x",i,sep="")}
  colnames(array) <- new_colnames
  if(time){colnames(array)[N+1] <- "time"} # Rename time column
  array # Return array
}

