
#=================================================================================#
#                       RATIO ANALYSIS OF EVOLUTION ARRAYS
#=================================================================================#


# Find the ratios between the steps for a given element array
element_ratios <- function(curr_element){
  # Get number of elements and steps
  M <- ncol(curr_element) - 2
  steps <- nrow(curr_element)
  # Initalize the stack
  ratio_stack <- rep(NA, M)
  for(i in 2:steps){
    # Get ratio of rows from current step 
    curr_ratios <- curr_element[i, 1:M]/curr_element[i-1, 1:M]
    # Stack
    ratio_stack <- rbind(ratio_stack, curr_ratios)
  }
  ratio_stack
}

# Append ratio of row elements by each step
append_ratios <- function(evolved_batch){
  # Extract B
  B <- max(evolved_batch$element_index)
  # Assuming two non-element columns (time, index) initialize the ratio stack
  r_stack <- element_ratios(element_array(evolved_batch, 1))
  for(i in 2:B){
    curr_ratios <- element_ratios(element_array(evolved_batch, i))
    r_stack <- rbind(r_stack, curr_ratios)
  }
  # Standardize the column names
  r_stack <- standardize_colnames(r_stack, prefix = "r_")
  # Return evolved batch with the ratios
  cbind(evolved_batch, r_stack)
}

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
