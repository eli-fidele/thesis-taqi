#=================================================================================#
#                               RATIO STEADINESS ANALYSIS
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

# Appends steadiness statistics to an evolved batch
.append_steadiness <- function(evolved_batch, P){
  # Extract the ratio array of the evolved batch
  ratios <- .extract_ratios(evolved_batch)
  ratios <- ratios[,2:ncol(ratios)]
  # Get the largest eigenvalue
  eigenvalues <- spectrum(P)
  lambda_1 <- .read_eigenvalue(eigenvalues, which.max(eigenvalues$Norm))
  steadiness <- do.call("rbind", lapply(X = 1:nrow(ratios), FUN = .steadiness_stats, ratios, lambda_1))
  colnames(steadiness) <- c("mean_stdns","max_stdns","lst_stdy","mst_stdy") # Rename columns
  # Once the steadiness analysis, cbind the statistics columns with the evolved batch and return
  cbind(evolved_batch[,1:(ncol(ratios)+2)], steadiness)
}

# Helper function
.steadiness_stats <- function(i, ratios, lambda_1){
  # Get current row of ratios
  curr <- ratios[i,1:ncol(ratios)]
  # If encountered a time = 0 row, return NAs
  if(is.na(curr[[1]])){
    rep(NA, 4)
  }
  else{  
    diff0 <- abs(curr - lambda_1) # Find absolute difference from largest eigenvalue
    diff1 <- abs(curr - abs(lambda_1)) # Find absolute difference from largest eigenvalue
    # Row statistics
    mean_stdns <- mean(as.numeric(diff1))
    max_stdns <- max(diff1)
    # Column statistics
    lst_stdy <- which.min(diff1)[[1]]
    mst_stdy <- which.max(diff1)[[1]]
    # Result
    c(mean_stdns, max_stdns, lst_stdy, mst_stdy)
  }
}