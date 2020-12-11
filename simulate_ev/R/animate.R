
#=================================================================================#
#                               PLOTTING FUNCTIONS
#=================================================================================#



#=================================================================================#
#                       EVOLUTION ARRAY ANALYSIS FUNCTIONS
#=================================================================================#

# Extract the array for a particular element
element_array <- function(evolved_batch, index){
  evolved_batch %>% filter(element_index == index)
}

# Extract the array for a range of elements
elements_array <- function(evolved_batch, indices){
  evolved_batch %>% filter(element_index %in% indices)
}

# Extract the array for a particular time
time_array <- function(evolved_batch, at_time){
  evolved_batch %>% filter(time == at_time)
}

# Extract the array for a range of time
times_array <- function(evolved_batch, time_range){
  evolved_batch %>% filter(time %in% time_range)
}

append_ratios_3d <- function(curr_array){
  r_x1 <- rep(1, nrow(curr_array))
  r_x2 <- rep(1, nrow(curr_array))
  r_x3 <- rep(1, nrow(curr_array))
  for(i in 2:nrow(curr_array)){
    # Find ratios between one step to the other for each variable
    r_x1[i] <- curr_array$x1[i] / curr_array$x1[i - 1]
    r_x2[i] <- curr_array$x2[i] / curr_array$x2[i - 1]
    r_x3[i] <- curr_array$x3[i] / curr_array$x3[i - 1]
  }
  cbind(curr_array, r_x1, r_x2, r_x3)
}

append_ratios_2d <- function(curr_array){
  r_x1 <- rep(1, nrow(curr_array))
  r_x2 <- rep(1, nrow(curr_array))
  for(i in 2:nrow(curr_array)){
    # Find ratios between one step to the other for each variable
    r_x1[i] <- curr_array$x1[i] / curr_array$x1[i - 1]
    r_x2[i] <- curr_array$x2[i] / curr_array$x2[i - 1]
  }
  cbind(curr_array, r_x1, r_x2)
}
  
# create wrangling functions that do the following
# removes initial element of the evolution array
# burns in every nth element (or keep in evol?)

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
evolve_batch <- function(batch, steps, burn_in = 1, with_steps = F){
  evol_stack <- evolve(batch[1,], P, steps, burn_in, with_steps) #append first batch element evolution array
  for(i in 2:B){ 
    evol <- evolve(batch[i,], P, steps, burn_in, with_steps) # obtain evol array of current row of the batch 
    evol_stack <- rbind(evol_stack, evol) # append rest of batch element arrays by stacking
  }
  rownames(evol_stack) <- 1:nrow(evol_stack) # standardize row names
  # Return the batch with added indices
  indexed_batch(data.frame(evol_stack), steps)
}

#=================================================================================#
#                               ELEMENTARY FUNCTIONS
#=================================================================================#

evolve <- function(v, P, steps, burn_in = 1, with_steps = F){
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
    labs(title = paste("Evolution of Monte Carlo Batch",mat_str))
}

#plots the evolution arrays of a 2d evolved batch
batch_2d_plot <- function(batch_data, mat_str = ""){
  ggplot(batch_data, mapping = aes(x = x1, y = x2, color = as.factor(element_index))) + 
    geom_point() +
    theme(legend.position = "none") +
    labs(title = paste("Evolution of 2D Monte Carlo Batch",mat_str))
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
  element_index <- rep(NA, nrow(evolved_batch))
  # append indexing by running loops
  for(i in 1:nrow(evolved_batch)){element_index[i] <- 1 + floor((i-1)/(steps+1))}
  # only append index column if it is not there already
  if(!("element_index" %in% colnames(evolved_batch))){evolved_batch <- cbind(evolved_batch,element_index)}
  evolved_batch
}