
#=================================================================================#
#                                 SIMULATION CODE
#=================================================================================#

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
  # Elements may be uniform only if the matrix isn't stochastic
  batch <- generate_batch(N = ncol(P), batch_size, stoch = .isStochastic(P))
  # Evolve the batch and return it
  evolved_batch <- evolve_batch(batch, P, steps)
  list(batch, evolved_batch)
}
