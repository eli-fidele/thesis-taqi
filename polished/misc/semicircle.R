
# Code to generate samples from Wigners Semicircle Distribution

# Sample n observations from the Wigner Semicircle Distribution of radius R
rsemicircle <- function(n = 1, R){
  # Simulate uniform theta from [0, pi] (angles spanning semicircle in complex plane)
  #thetas <- pi * runif(n, min = 0, max = 1)  
  # Normalize units to be of hypotenuse length at theta on a circle of radius R
  #sample <- (R * cos(pi/2 - thetas))
  #sample <- (R * cos(pi - thetas))
  x_comps <- runif(n, -R, R)
  sample <- (sqrt(R^2 - x_comps^2))
}

# Sample instance of the semicircle distribution
# r_semicircle<-function(R){
#   max_y <- 2/(pi*R) # Maximum of density function
#   repeat {
#     x <- runif(1, 0, 1)
#     y <- runif(1, 0, max_y)
#     fx <- (2/(pi * R^2)) * sqrt(R^2 - x^2) # Semicircle density function
#     if (y < fx) 
#       return(x)
#   }
# }

# rz2<-function (n) 
# {
#   zvector <- vector("numeric", n)
#   for (i in 1:n) {
#     zvector[i] <- r_semicircle()
#   }
#   zvector
# }

r_semicircle <- function(n, R){
  max_y <- 2/(pi*R) # Maximum of density function
  # Sample xs and ys
  xvec <- runif(n, -R, R)
  yvec <- runif(n, 0, max_y)
  # Get densities and reject
  fvec <- (2/(pi * R^2)) * sqrt(R^2 - xvec^2) # Semicircle density function
  xvec[yvec < fvec]
}

sample_semicircle <- function(n, R){
  x <- r_semicircle(n, R)
  len <- length(x)
  aprob <- len/n
  shortby <- n - len
  n2 <- round(shortby/aprob)
  x2 <- r_semicircle(n2, R)
  x3 <- c(x, x2)
  #c(-x3,x3)
  x3
}