
# This script includes functions that are metrics for analyzing eigenvectors.
# The inputs for these metric functions have been standardized to the dataframe object returned by the eigen_frame() function
# in eigenvectors.R.

#=================================================================================#
#                               EIGENVECTOR METRICS
#=================================================================================#

# Returns a tidied dataframe of the eigenvalues of a random matrix
eigen_frame <- function(P){
  M <- length(P[1,])
  eigenvectors <- data.frame(eigen(P)[2])
  complex <- matrix(rep(NA,3*M*M), ncol = 3) # set 3 to hold (re,im) pair and whose row it belongs to
  colnames(complex) <- c("Re","Im","row_i")
  for(i in 1:M){
    for(j in 1:M){
      curr <- eigenvectors[i,j]
      complex[ M*(i-1) + j, ] <- c(round(Re(curr),5),round(Im(curr),5),i) 
    }
  }
  data.frame(complex)
}
