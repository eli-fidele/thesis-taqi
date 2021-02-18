

#=================================================================================#
#                             ENSEMBLE ANALYSIS
#=================================================================================#

# Combines the results of a ensemble simulation into a master array
glue_arrays <- function(ensemble_sim, array_index = 1){
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
  
#=================================================================================#
#                             RANDOM MATRIX ENSEMBLES 
#=================================================================================#

# Random matrix ensemble
## FIX: generalize to multiple arguments rather than just dimensions, accounting for defaults
RM_ensemble <- function(mat_type, args, size){
  # If the mat_type entry is not the function itself, parse the mat_type string
  if(class(mat_type) == "function"){
    fxn <- mat_type
  } else{
    fxn <- parse_RMfxn(mat_type)  
  }
  # Replicate with appopriate number of arguments by using parse_args
  ensemble <- replicate(n = size, expr = parse_args(fxn, args), simplify = F)
  # Return the ensemble
  ensemble
}

#=================================================================================#
#                                  LOW-LEVEL 
#=================================================================================#


parse_args <- function(fxn, args){
  # Get number of arguments
  n_args <- length(args)
  if(n_args == 1){return(fxn(args[[1]]))}
  if(n_args == 2){return(fxn(args[[1]], args[[2]]))}
  if(n_args == 3){return(fxn(args[[1]], args[[2]], args[[3]]))}
  if(n_args == 4){return(fxn(args[[1]], args[[2]], args[[3]], args[[4]]))}
  if(n_args == 5){return(fxn(args[[1]], args[[2]], args[[3]], args[[4]], args[[5]]))}
  if(n_args == 6){return(fxn(args[[1]], args[[2]], args[[3]], args[[4]], args[[5]], args[[6]]))}
}

# Take a list of arguments, and an index holding a vector to properly return a vector of the arguments
# Steps: get the length of the vector argument, then do args[[i]][j] for j in 1:len(vec)
vector_arg <- function(args, index){
  NA
}

# Elementary matrix type string parser
# FUTURE: develop to enable string subsets to be appopriately parsed.
parse_RMfxn <- function(type_str){
  # Basic lexicon of string interpretation
  normal_str <- c("Normal", "normal","norm","n","N")
  stoch_str <- c("Stochastic","stochastic", "stoch", "st", "s", "S")
  erdos_str <- c("Erdos","erdos","er","e", "ER", "E")
  # Test whether its normal
  if(type_str %in% normal_str){
    fxn <- RM_normal
  }
  # Test whether its stochastic
  if(type_str %in% stoch_str){
    fxn <- RM_stoch
  }
  # Test whether its an ER-graph stochastic
  if(type_str %in% erdos_str){
    fxn <- RM_erdos
  }
  # If all conditions fail, return NA
  else{
    fxn <- NA
    print("Please try again.")
  }
  # Return function
  fxn
}

#=================================================================================#
#                                  EXPERIMENTAL 
#=================================================================================#

# Arbitrary stacker function which returns a stacked array of the results
# Dimensions to ensure correct sizing, or unnecessary?
# Argument array conventions?
stacker <- function(fxn, args, iter, dimensions){
  stack <- fxn(args[1,]) # Initialize the array
  for(i in 2:iter){ # Iterate over the rest of the argument array
    curr <- fxn(args[i,])
    stack <- rbind(stack, curr)
  }
  stack # Return the stack
}

