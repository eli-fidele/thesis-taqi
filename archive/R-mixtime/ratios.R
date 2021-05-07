
#=================================================================================#
#                       EIGENVALUE CLASSIFICATION OF RATIOS
#=================================================================================#

# Takes in an **FINAL-TIME** evolved batch to return an analysis of the simluated eigenvalues
.eigen_classify <- function(evolved_batch, P, epsilon = 0.1){
  # Extract the ratio array of the evolved batch
  ratios <- .extract_ratios(evolved_batch)
  # Get the eigenvalue array
  eigenvalues <- eigen(P)$values
  # Create associated eigenvalue index column (with 0 implying no match)
  eigen_index <- rep(0, nrow(evolved_batch))
  # Sift through the numerical "eigenvalues" and find potential matches
  for(i in 1:nrow(ratios)){
    # Go row by row, ignoring the element index column
    curr <- ratios[i,2:ncol(ratios)] 
    # Try to see if any of the eigenvalues fit
    for(K in 1:length(eigenvalues)){
      curr_eigenvalue <- eigenvalues[K] # Load current eigenvalue as a numerical type for comparison (function in eigen.R)
      if(.candidate_eigenvector(curr, curr_eigenvalue, epsilon)){
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
.candidate_eigenvector <- function(vector, eigenvalue, epsilon){
  # Componentwise check on whether each ratio is within epsilon of the eigenvalue
  boolean_vector <- (abs(vector - eigenvalue) < epsilon)
  # If all the values are true, the vector is a candidate
  !((FALSE %in% boolean_vector) || (NA %in% boolean_vector))
}

#=================================================================================#
#                         EVOLUTION ARRAY RATIO-FILTERING
#=================================================================================#

# Helper function, returns a vector of all the ratio entries at a given time
ratios_by_time <- function(evolved_batch, at_time, log = T){
  # Extract ratios for a given time and remove the 'element_index' column.
  ratios <- .extract_ratios(by.time(evolved_batch, at_time))
  ratios <- ratios[,2:ncol(ratios)] # Drop non-ratio columns
  if(log){ratios <- log(ratios)} # Take log if prompted
  all_ratios <- as.vector(ratios$r_x1) # Initialize vector by taking ratios in first row
  for(col in 2:ncol(ratios)){ 
    curr_row <- as.vector(ratios[,col]) 
    all_ratios <- c(all_ratios, curr_row) # Concatenate the rest of the consective rows' ratios
  }
  all_ratios # Return ratios
}

# Extracts the ratio columns for a given batch
.extract_ratios <- function(evolved_batch){
  # Get number of dimensions and steps (two for time and element_index; divide by two for xi and r_xi)
  M <- (ncol(evolved_batch) - 2)/2
  # Return the array with only the ratios, assumming two non-element columns        
  evolved_batch[,(M+2):(2*M+2)]
}

#=================================================================================#
#                    RATIO VARIANCE ANALYSIS OF EVOLUTION ARRAYS
#=================================================================================#

# Gives the variance of the ratio entries for all the columns by time
variance_by_time <- function(evolved_batch, at_time, log = T){
  variances <- rep(NA, length(at_time)) # Create a vector to hold the variance for each time
  for(i in 1:length(at_time)){ 
    curr_var <- var(ratios_by_time(evolved_batch, at_time = i, log), na.rm = T) # Get the variance at that time
    #if(log){curr_var <- log(curr_var)} # Take log if prompted
    variances[i] <- curr_var
  }
  variances # Return variances
}

