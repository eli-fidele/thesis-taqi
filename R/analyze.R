
#=================================================================================#
#                       EVOLUTION ARRAY FILTERING FUNCTIONS
#=================================================================================#

# Takes in an evolved batch to return an analysis of the eigenvalue mixing times
eigen_mixtime <- function(evolved_batch, P, epsilon = 0.1){
  # Extract the ratio array of the evolved batch
  ratios <- extract_ratios(evolved_batch)
  # Get the eigenvalue array
  eigenvalues <- eval_frame(P)
  # Create associated eigenvalue index column (with 0 implying no match)
  eigen_index <- rep(0, nrow(evolved_batch))
  # Sift through the numerical "eigenvalues" and find potential matches
  for(i in 1:nrow(ratios)){
    # Go row by row, ignoring the element index column
    curr <- ratios[i,2:ncol(ratios)] 
    # Try to see if any of the eigenvalues fit
    for(K in 1:nrow(eigenvalues)){
      curr_eigenvalue <- read_eigenvalue(eigenvalues, K) # Load current eigenvalue as a numerical type for comparison (function in eigen.R)
      if(candidate_eigenvector(curr, curr_eigenvalue, epsilon)){
        # If a candidate eigenvector is found, record this in the index column and exit the loop
        eigen_index[i] <- K
        break
      }
    }
  }
  # Once the eigenvalue crosscheck is complete, cbind the Index column with the evolved batch
  cbind(evolved_batch, eigen_index)
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
