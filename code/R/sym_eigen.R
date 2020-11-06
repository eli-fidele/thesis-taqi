
#=================================================================================#
#                           RANDOM SYMMETRIC MATRICES
#=================================================================================#

RM_symm <- function(M,f,ep){
  dist <- data.frame(x = runif(M**2, ep*(f-1), ep*f))
  make_symm(dist)
}

#=================================================================================#
#                     RANDOM SYMMETRIC MATRICES HELPER FXNS
#=================================================================================#

#returns a vector with M^2 entries that are uniformly distributed with a fraction of its values positive = f 
unif_frand <- function(M,f=T,ep){
  # unless specifically initialized, a random fraction will be chosen
  if(f){
    f <- runif(1,0,1)
    paste("f: ",f,sep="")
  }
  dist <- data.frame(x = runif(M**2, ep*(f-1), ep*f))
  dist <- dist %>% mutate(x_pos = ifelse(x > 0, 1, 0))
  dist
}

make_symm <- function(dist){
  N <- sqrt(length(dist$x))
  P <- matrix(data = dist$x, nrow = N, ncol = N)
  P[lower.tri(P)] <- P[upper.tri(P)]
  P
}

# returns proportion of positive entries of any matrix P
pos_entries <- function(P){
  pos_entries <- length(matrix(P[P[,] > 0], nrow = 1))
  pos_entries/(length(P))   
}