
#=================================================================================#
#                             ENSEMBLE SIMULATION 
#=================================================================================#

ens_mixtime <- function(ensemble, batch_size, steps, epsilon = 0.1){
  # Get a sample matrix
  P <- ensemble[[1]]
  # Batch elements may be uniform only if the matrix isn't stochastic
  initial <- generate_batch(N = ncol(P), batch_size = batch_size, stoch = .isStochastic(P))
  # Evolve every matrix in the ensemble
  #do.call("c",lapply(ensemble, get_mixtime, initial, steps, epsilon = 0.1))
  map_dfr(ensemble, get_mixtime, initial, steps, epsilon = 0.1)
}

#=================================================================================#
#                       MIXTIME SIMULATION WRAPPER FUNCTIONS
#=================================================================================#

get_mixtime <- function(P, initial, steps, epsilon = 0.1, simple = T){
  # Evolve the batch
  evolved <- evolve_batch(batch = initial, P, steps)
  # Add the eigen index to the evolved batch
  classified <- .eigen_classify(evolved, P, epsilon) # Classify candidate eigenvectors
  # Using the eigen index, classify the initial elements and get their mixing times 
  eigen_index <- by.time(classified, at_time = steps)$eigen_index
  initial <- cbind(initial, eigen_index) 
  initial <- .eigen_mixtime(classified, initial, complete = FALSE)
  if(simple) { return(data.frame(mixtime = initial)) }
  # Return a list of the two arrays
  list(initial, evolved)
}
  

#=================================================================================#
#                       MIXTIME SIMULATION WRAPPER FUNCTIONS
#=================================================================================#

# Simulation to determine the eigenvalue mixing time of a random matrix.
sim_mixtime <- function(P, batch_size, steps, epsilon = 0.1){
  B <- batch_size
  init_sim <- sim_initial(P, B, steps) # Initialize the simulation
  # Extract the arrays
  initial <- init_sim[[1]]
  evolved <- init_sim[[2]]
  # Return list of arrays
  get_mixtime(initial, evolved, steps, epsilon)
}

# Generate and evolve a batch of points for a given random matrix P. 
# This function is a basic "initial" simulation. Other simulation functions will utilize this function as a base.
sim_initial <- function(P, batch_size, steps){
  # Batch elements may be uniform only if the matrix isn't stochastic
  initial <- generate_batch(N = ncol(P), batch_size = batch_size, stoch = .isStochastic(P))
  # Evolve the batch and return it
  evolved <- evolve_batch(initial, P, steps)
  list(initial, evolved)
}


