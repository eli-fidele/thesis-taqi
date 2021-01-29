
#=================================================================================#
#                            BATCH GENERATION FUNCTIONS
#=================================================================================#

# Run a batch
run_batch <- function(P, B = 50, lambda = 1, steps = 25, indexed = T, ratios = T){
  M <- ncol(P)
  # Make the batch
  batch <- make_batch(M, B)
  # Evolve the batch
  evolved_batch <- evolve_batch(batch, steps, lambda, with_steps = T)
  # Add indexing if prompted
  if(indexed){
    evolved_batch <- indexed_batch(evolved_batch, steps)
    # Since ratios are only possible with indexing, we nest the ratios argument here
    if(ratios){
      evolved_batch <- append_ratios(evolved_batch)}
    }
  evolved_batch
}

# Evolve each element of the batch vector a given number of steps 
evolve_batch <- function(batch, steps, burn_in = 1, with_steps = T){
  evol_stack <- evolve(batch[1,], P, steps, burn_in, with_steps) #append first batch element evolution array
  B <- nrow(batch)
  for(i in 2:B){ 
    evol <- evolve(batch[i,], P, steps, burn_in, with_steps) # Obtain evolution array of current element of the batch 
    evol_stack <- rbind(evol_stack, evol) # Append the rest of batch element arrays by stacking
  }
  rownames(evol_stack) <- 1:nrow(evol_stack) # standardize row names
  data.frame(evol_stack)
}

# Generate a Monte Carlo batch
make_batch <- function(M, B, lambda = 1){
  batch <- matrix(rep(NA, B * M), nrow = B)  # create [B x M] batch matrix
  # Generate the batch 
  for(i in 1:B){batch[i,] <- runif(n = M, min = -lambda, max = lambda)}
  batch <- standardize_colnames(batch) # standardize the column names
  data.frame(batch) # return batch
}

# Evolves an element of a batch by the matrix P
evolve <- function(v, P, steps, burn_in = 1, with_steps = T){
  M <- ncol(P)
  # Simulate the evolution matrix of a given batch element
  vals <- matrix(rep(NA, M * steps), ncol = M)
  vals <- standardize_colnames(vals)
  if(with_steps){vals <- cbind(vals, rep(0, steps+1))}
  # Evolve the batch element
  for(i in 1:steps){
    # Evolved vector (v * P^steps)
    evolved_row <- as.numeric(v) %*% matrix.power(P,burn_in*i)
    # Add time index to evolved row
    vals[i, ] <- cbind(evolved_row, i)
  }
  #store the values in a dataframe
  if(with_steps){
    vals <- rbind(c(as.numeric(v),0),vals)
  } else{
    vals <- rbind(as.numeric(v),vals)}
  # Rename steps column to 'time'
  colnames(vals)[ncol(vals)] <- "time"
  vals
}


#=================================================================================#
#                            INDEXING/NAMING FUNCTIONS
#=================================================================================#

standardize_colnames <- function(array, prefix = ""){
  # get string vector
  str_vec <- rep(NA, ncol(array))
  # rename the columns
  for(i in 1:ncol(array)){
    str_vec[i] <- paste(prefix,"x",i,sep="")
  }
  colnames(array) <- str_vec
  # return renamed col array
  array
}

# This method adds an index to clarify which rows belong to which batch element's evolution array it belongs
indexed_batch <- function(evolved_batch, steps){
  # create index column
  element_index <- rep(NA, nrow(evolved_batch))
  # append indexing by running loops
  for(i in 1:nrow(evolved_batch)){element_index[i] <- 1 + floor((i-1)/(steps+1))}
  # only append index column if it is not there already
  if(!("element_index" %in% colnames(evolved_batch))){evolved_batch <- cbind(evolved_batch,element_index)}
  evolved_batch
}

#=================================================================================#
#                               PLOTTING FUNCTIONS
#=================================================================================#

#plots the evolution arrays of a 2d evolved batch
batch_3d_plot <- function(batch_data,mat_str=""){
  plot_empty <- ggplot() 
  plot_12 <- batch_2d_customplot(batch_data, 1, 2,mat_str)
  plot_23 <- batch_2d_customplot(batch_data, 2, 3,mat_str)
  plot_13 <- batch_2d_customplot(batch_data, 1, 3,mat_str)
  grid.arrange(plot_empty,plot_12,plot_23,plot_13, ncol = 2)
}

#plots the evolution arrays of a 2d evolved batch
batch_2d_customplot <- function(batch_data, n1, n2, mat_str = ""){
  ggplot(batch_data, mapping = aes(color = as.factor(element_index))) + 
    geom_point(mapping = aes_string(x = paste("x",n1,sep=""), y = paste("x",n2,sep=""))) +
    theme(legend.position = "none") +
    labs(title = paste("Evolution of a Markov Chain",mat_str))
}

#plots the evolution arrays of a 2d evolved batch
batch_2d_plot <- function(batch_data, mat_str = ""){
  ggplot(batch_data, mapping = aes(x = x1, y = x2, color = as.factor(element_index))) + 
    geom_point() +
    theme(legend.position = "none") +
    labs(title = paste("Evolution of a Markov Chain",mat_str))
}