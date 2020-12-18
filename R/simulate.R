
#=================================================================================#
#                               PLOTTING FUNCTIONS
#=================================================================================#

extract_evol_array <- function(evolved_batch, index){
  evolved_batch %>% filter(index_column == index)
}

#=================================================================================#
#                       EVOLUTION ARRAY ANALYSIS FUNCTIONS
#=================================================================================#

# create wrangling functions that do the following
# removes initial element of the evolution array
# burns in every nth element (or keep in evol?)

#=================================================================================#
#                            BATCH GENERATION FUNCTIONS
#=================================================================================#

# Run a batch
run_batch <- function(P, B = 50, lambda = 1, steps = 25){
  M <- ncol(P)
  batch <- make_batch(M, B)
  evolve_batch(batch, steps, lambda)
}

# Evolve each element of the batch vector a given number of steps 
evolve_batch <- function(batch, steps, burn_in = 1, with_steps = T){
  evol_stack <- evolve(batch[1,], P, steps, burn_in, with_steps) #append first batch element evolution array
  B <- nrow(batch)
  for(i in 2:B){ 
    evol <- evolve(batch[i,], P, steps, burn_in, with_steps) # obtain evol array of current row of the batch 
    evol_stack <- rbind(evol_stack, evol) # append rest of batch element arrays by stacking
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

#=================================================================================#
#                               ELEMENTARY FUNCTIONS
#=================================================================================#

evolve <- function(v, P, steps, burn_in = 1, with_steps = F){
  it <- steps
  M <- ncol(P)
  # simulate and record evolution of pi
  vals <- matrix(rep(NA, M * it), ncol = M)
  vals <- standardize_colnames(vals)
  # evolve pi 
  for(i in 1:it){
    vals[i, ] <- as.numeric(v) %*% matrix.power(P,burn_in*i)
  }
  #store the values in a dataframe
  vals <- rbind(v,vals)
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
  ggplot(batch_data, mapping = aes(color = as.factor(index_column))) + 
    geom_point(mapping = aes_string(x = paste("x",n1,sep=""), y = paste("x",n2,sep=""))) +
    theme(legend.position = "none") +
    labs(title = paste("Evolution of Monte Carlo Batch",mat_str))
}

#plots the evolution arrays of a 2d evolved batch
batch_2d_plot <- function(batch_data, mat_str = ""){
  ggplot(batch_data, mapping = aes(x = x1, y = x2, color = as.factor(index_column))) + 
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
  index_column <- rep(NA, nrow(evolved_batch))
  # append indexing by running loops
  for(i in 1:nrow(evolved_batch)){index_column[i] <- 1 + floor((i-1)/(steps+1))}
  # only append index column if it is not there already
  if(!("index_column" %in% colnames(evolved_batch))){evolved_batch <- cbind(evolved_batch,index_column)}
  evolved_batch
}