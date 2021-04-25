
#=================================================================================#
#                                 SIMULATION CODE
#=================================================================================#

# Simulation for analyzing the steadiness statistics of the consecutive ratio sequence
sim_steadiness <- function(P, batch_size, steps){
  init_sim <- sim_initial(P, batch_size, steps) # Initialize the simulation
  # Extract the arrays
  batch <- init_sim[[1]]
  evolved_batch <- init_sim[[2]]
  # Perform the steadiness analysis
  evolved_batch <- .append_steadiness(evolved_batch, P)
  # Return list of arrays
  list(batch, evolved_batch)
}

# Simulation to determine the eigenvalue mixing time of a random matrix.
sim_mixtime <- function(P, batch_size, steps, epsilon = 0.1){
  B <- batch_size
  init_sim <- sim_initial(P, B, steps) # Initialize the simulation
  # Extract the arrays
  batch <- init_sim[[1]]
  evolved_batch <- init_sim[[2]]
  # Perform the mixing time analysis
  evolved_batch <- .eigen_classify(evolved_batch, P, epsilon) # Classify candidate eigenvectors
  eigen_index <- by.time(evolved_batch, at_time = steps)$eigen_index # Get the array at the final time
  batch <- cbind(batch, eigen_index) # Add the eigen_index of the batch elements to the array
  batch <- .eigen_mixtime(evolved_batch, batch) # Add the mixing times to the batch element array
  # Return list of arrays
  list(batch, evolved_batch)
}

# Generate and evolve a batch of points for a given random matrix P. 
# This function is a basic "initial" simulation. Other simulation functions will utilize this function as a base.
sim_initial <- function(P, batch_size, steps){
  # Batch elements may be uniform only if the matrix isn't stochastic
  batch <- generate_batch(N = ncol(P), batch_size = batch_size, stoch = .isStochastic(P))
  # Evolve the batch and return it
  evolved_batch <- evolve_batch(batch, P, steps)
  list(batch, evolved_batch)
}

#=================================================================================#
#                             ENSEMBLE ANALYSIS
#=================================================================================#

# Combines the results of a ensemble simulation into a master array
sim.glue_arrays <- function(ensemble_sim, array_index = 1){
  batch <- ensemble_sim[[1]][[array_index]] # Get first result
  batch_size <- nrow(batch) # Get batch size for reference
  mat_idx <- data.frame(mat_idx = rep(1, batch_size)) # Create matrix index column
  batch <- cbind(batch, mat_idx) # Initialize the master batch
  # Repeat for the rest of the elements
  for(i in 2:length(ensemble_sim)){
    curr_batch <- ensemble_sim[[i]][[array_index]] # Get array for matrix i 
    mat_idx <- data.frame(mat_idx = rep(i, batch_size)) # Index matrix i
    batch <- rbind(batch, cbind(curr_batch, mat_idx)) # Concatenate arrays
  }
  batch
}

#=================================================================================#
#                             ENSEMBLE SIMULATION 
#=================================================================================#

# Simulates the mixtimes for an ensemble of matrices
mixtime_ensemble <- function(ensemble, batch_size, steps, epsilon = 0.1){
  # Initialize the stack
  ensemble_result <- list(sim_by_element(ensemble, batch_size, steps, epsilon, ensemble_index = 1))
  # Go through rest of ensemble
  for(i in 2:length(ensemble)){
    curr_result <- sim_by_element(ensemble, batch_size, steps, epsilon, ensemble_index = i)
    ensemble_result <- c(ensemble_result, list(curr_result)) # Concatenate results
  }
  ensemble_result
}

sim_by_element <- function(ensemble, batch_size, steps, epsilon, ensemble_index){
  P <- ensemble[[ensemble_index]] # Extract the matrix
  sim <- mixtime_sim(P, batch_size, steps, epsilon) # Get the simulation list for one matrix
  c(sim, list(P))  # Extract the results alongside the matrix
}