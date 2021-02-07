

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
#                               BATCH EXECUTION 
#=================================================================================#


# Executes a batch of run_batch function on elements of a random matrix ensemble.
execute_batch <- function(RM_mattype, RM_args, batch_args = list(100), execution_class, stack_args = list(20, TRUE)){
  # Create the random matrix ensemble
  ensemble <- RM_ensemble(RM_mattype, RM_args, 10)
  # Obtain variables
  M <- ncol(ensemble[[1]])
  # Create the batch
  batch <- make_batch(M, B = batch_args[[1]])
  # Stack the batches
  stack <- execution_class(ensemble[[1]], batch)
  for(i in 2:length(ensemble)){
    curr <- execution_class(ensemble[[i]], batch, stack_args)
    stack <- rbind(stack, curr)
  }
  # Return the evolved stack
  stack
}

run_class0 <- function(P, batch, stack_args = list(20, TRUE)){
  # Unwind the variables
  steps <- stack_args[[1]]
  with_ratios <- stack_args[[2]]
  # Evolve the batch and return it
  evolve_batch(P, batch, steps, with_ratios)
}

#=================================================================================#
#                                  LOW-LEVEL 
#=================================================================================#


parse_args <- function(fxn, args){
  # Get number of arguments
  n_args <- length(args)
  if(n_args == 1){fxn(args[[1]])}
  if(n_args == 2){fxn(args[[1]], args[[2]])}
  if(n_args == 3){fxn(args[[1]], args[[2]], args[[3]])}
  if(n_args == 4){fxn(args[[1]], args[[2]], args[[3]], args[[4]])}
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

