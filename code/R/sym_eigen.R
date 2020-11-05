
#=================================================================================#
#                           RANDOM SYMMETRIC MATRICES
#=================================================================================#

RM_symm <- function(M,f,ep){
  M <- 3
  f <- 0.5
  ep <- 100
  dist <- unif_fpos(M = M, f = f, ep = ep)
  make_symm(dist)
}

r_normal <- function(M, mu, sd){
  rnorm(n = M, mean = mu, sd = sd)
}

#=================================================================================#
# 
#=================================================================================#

#returns a vector with M^2 entries that are uniformly distributed with a fraction of its values positive = f 
unif_fpos <- function(M,f,ep){
  # unless specifically initialized, a random fraction will be chosen
  if(F){
    f <- runif(1,0,1)
    paste("f: ",f,sep="")
  }
  b <- f
  a <- (f-1)
  dist <- data.frame(x = runif(M**2, ep*a, ep*b))
  dist <- dist %>% mutate(x_neg = ifelse(x < 0,yes = 1, no = 0))
  dist
}

make_symm <- function(dist){
  N <- sqrt(length(dist$x))
  P <- matrix(data = dist$x, nrow = N, ncol = N)
  LT <- lower.tri(P)
  UT <- upper.tri(P)
  P[LT] <- P[UT]
  P
}

# returns proportion of positive entries of any matrix P
pos_entries <- function(P){
  pos_entries <- length(matrix(P[P[,] > 0], nrow = 1))
  pos_entries/(length(P))   
}