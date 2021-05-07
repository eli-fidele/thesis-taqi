#=================================================================================#
#                               MIXING TIME ANALYSIS
#=================================================================================#

.eigen_mixtime <- function(classified_batch, batch, complete = F){
  # Extract B (number of elements)
  B <- max(classified_batch$element_index)
  # Create initial mixtime vector
  mixtime <- rep(0, B)
  # Loop over every element of the batch, finding the mixing time
  for(i in 1:B){
    # Extract the classified column of eigen_indices over time
    seq <- by.element(classified_batch, index = i)$eigen_index
    # Find the time such that the eigen_index is non-zero (implying near-convergence)
    mixtime[i] <- min(which(seq != 0)) - 1
  }
  mixtime[which(mixtime == Inf)] <- NA # Address Inf entries by NA'ing them
  if(complete){ return (cbind(batch, mixtime)) }
  # Return the mixing times
  mixtime  
}

# Outputs the proportion of batch elements which are mixed
.prop_mixed <- function(batch){
  prop_unmixed <- 1 - length(which(is.na(batch$mixtime)))/nrow(batch)
  paste("This batch is ",100*round(prop_unmixed, 3),"% mixed.",sep="")
}