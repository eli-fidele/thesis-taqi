
#=================================================================================#
#                             BATCH ANALYSIS FUNCTIONS
#=================================================================================#



#=================================================================================#
#                            BATCH GENERATION FUNCTIONS
#=================================================================================#

# Generate a Monte Carlo batch
make_batch <- function(M, B, lambda = 1){
  batch <- matrix(rep(NA, B * M), nrow = B)  # create [B x M] batch matrix
  # Generate the batch 
  for(i in 1:B){batch[i,] <- runif(n = M, min = -lambda, max = lambda)}
  batch <- standardize_colnames(batch) # standardize the column names
  data.frame(batch) # return batch
}

# Evolve each element of the batch vector a given number of steps 
evolve_batch <- function(batch, steps){
  evol_stack <- evolve(batch[1,], P, steps) #append first batch element evolution array
  for(i in 2:B){ 
    evol <- evolve(batch[i,], P, steps) # obtain evol array of current row of the batch 
    evol_stack <- rbind(evol_stack, evol) # append rest of batch element arrays by stacking
  }
  rownames(evol_stack) <- 1:nrow(evol_stack) # standardize row names
  data.frame(evol_stack)
}

#=================================================================================#
#                               ELEMENTARY FUNCTIONS
#=================================================================================#

evolve <- function(v, P, steps){
  it <- steps
  M <- ncol(P)
  # simulate and record evolution of pi
  vals <- matrix(rep(NA, M * it), ncol = M)
  vals <- standardize_colnames(vals)
  # evolve pi 
  for(i in 1:it){
    vals[i, ] <- as.numeric(v) %*% matrix.power(P,i)
  }
  #store the values in a dataframe
  rbind(v,vals)
}

distance <- function(pi,ref_dist){
  #plot difference from a reference/stationary distribution
  diff <- rbind(evolve(pi),ref_dist)
  dist_vec <- rep(0, it)
  for(i in 1:it){
    curr_dist <- stats::dist(diff[c(i,it+1),], method = "euclidean")
    dist_vec[i] <- curr_dist
  }
  data.frame(dist_vec)
}

#=================================================================================#
#                               HELPER FUNCTIONS
#=================================================================================#

standardize_colnames <- function(array){
  # get string vector
  str_vec <- rep(NA, ncol(array))
  # rename the columns
  for(i in 1:ncol(array)){
    str_vec[i] <- paste("x",i,sep="")
  }
  colnames(array) <- str_vec
  # return renamed col array
  array
}

# This method adds an index to clarify which rows belong to which batch element's evolution array it belongs
indexed_batch <- function(evolved_batch, steps){
  # create index column
  index_column <- rep(NA, nrow(evolved_batch))
  # append indexing by running loops
  for(i in 1:nrow(evolved_batch)){index_column[i] <- 1 + floor((i-1)/(steps+1))}
  cbind(evolved_batch,index_column)
}