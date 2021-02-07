
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
    seq <- by_element(evolved_batch, index = i)$eigen_index
    # Find the time such that the eigen_index is non-zero (implying near-convergence)
    mixtime[i] <- min(which(seq != 0)) - 1
  }
  mixtime[which(mixtime == Inf)] <- NA # Address Inf entries by NA'ing them
  cbind(batch, mixtime)
}

#=================================================================================#
#                       EIGENVALUE ANALYSIS OF EVOLUTION ARRAYS
#=================================================================================#

# Takes in an **FINAL-TIME** evolved batch to return an analysis of the simluated eigenvalues
eigen_classify <- function(evolved_batch, P, epsilon = 0.1){
  # Extract the ratio array of the evolved batch
  ratios <- extract_ratios(evolved_batch)
  # Get the eigenvalue array
  eigenvalues <- spectrum(P)
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

# Given a ratio row, check to see if a vector is a candidate eigenvector
candidate_eigenvector <- function(vector, eigenvalue, epsilon){
  # Componentwise check on whether each ratio is within epsilon of the eigenvalue
  boolean_vector <- (abs(vector - eigenvalue) < epsilon)
  # If all the values are true, the vector is a candidate
  !((FALSE %in% boolean_vector) || (NA %in% boolean_vector))
}

#=================================================================================#
#                       RATIO ANALYSIS OF EVOLUTION ARRAYS
#=================================================================================#


# Find the ratios between the steps for a given element array
element_ratios <- function(curr_element){
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

# Append ratio of row elements by each step
append_ratios <- function(evolved_batch){
  # Extract B
  B <- max(evolved_batch$element_index)
  # Assuming two non-element columns (time, index) initialize the ratio stack
  r_stack <- element_ratios(by_element(evolved_batch, 1))
  for(i in 2:B){
    curr_ratios <- element_ratios(by_element(evolved_batch, i))
    r_stack <- rbind(r_stack, curr_ratios)
  }
  # Standardize the column names
  r_stack <- standardize_colnames(r_stack, prefix = "r_")
  # Return evolved batch with the ratios
  cbind(evolved_batch, r_stack)
}

#=================================================================================#
#                       EVOLUTION ARRAY FILTERING FUNCTIONS
#=================================================================================#

# Extracts the ratio columns for a given batch
extract_ratios <- function(evolved_batch){
  # Get number of dimensions and steps (two for time and element_index; divide by two for xi and r_xi)
  M <- (ncol(evolved_batch) - 2)/2
  # Return the array with only the ratios
  evolved_batch[,(M+2):(2*M+2)]
}

# Extract the array for a particular element/a range of elements
by_element <- function(evolved_batch, index){
  if(class(index) == "numeric"){
    evolved_batch %>% filter(element_index == index)
  } else {
    evolved_batch %>% filter(element_index %in% index)
  }
}

# Extract the array for a particular time/a range of times
by_time <- function(evolved_batch, at_time){
  if(class(at_time) == "numeric"){
    evolved_batch %>% filter(time == at_time)
  } else {
    evolved_batch %>% filter(time %in% at_time)
  }
}

