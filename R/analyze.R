
#=================================================================================#
#                               MIXING TIME ANALYSIS
#=================================================================================#

eigen_mixtime <- function(evolved_batch, batch){
  # Extract B (number of elements)
  B <- max(evolved_batch$element_index)
  # Create initial mixtime vector
  mixtime <- rep(0, B)
  # Loop over every element of the batch, finding the mixing time
  for(i in 1:B){
    # Extract the classified column of eigen_indices over time
    seq <- element_array(evolved_batch, index = i)$eigen_index
    # Find the time such that the eigen_index is non-zero (implying near-convergence)
    mixtime[i] <- min(which(seq != 0)) - 1
  }
  mixtime[which(mixtime == Inf)] <- NA # Address Inf entries by NA'ing them
  cbind(batch, mixtime)
}

#=================================================================================#
#                       EVOLUTION ARRAY FILTERING FUNCTIONS
#=================================================================================#

# Extracts the ratio columns for a given batch
extract_ratios <- function(evolved_batch){
  # Get number of dimensions and steps (two for time and element_index; divide by two for xi and r_xi)
  M <- (ncol(evolved_batch) - 2)/2
  # Return the array with only the ratios
  evolved_batch[,(M+2):(2*M + 2)]
}

# Extract the array for a particular element
element_array <- function(evolved_batch, index){
  evolved_batch %>% filter(element_index == index)
}

# Extract the array for a range of elements
elements_array <- function(evolved_batch, indices){
  evolved_batch %>% filter(element_index %in% indices)
}

# Extract the array for a particular time
time_array <- function(evolved_batch, at_time){
  evolved_batch %>% filter(time == at_time)
}

# Extract the array for a range of time
times_array <- function(evolved_batch, time_range){
  evolved_batch %>% filter(time %in% time_range)
}
