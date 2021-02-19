
#=================================================================================#
#                           SIMULATION FUNCTIONS
#=================================================================================#

#### SYMM_EVC_M (BEST FILE)

simulate_by_f_PROP <- function(M_max,f,ep,draws){
  M_vec <- sample(1:M_max, draws, replace = F)
  table <- data.frame(M = M_vec)
  
  prop_vec <- rep(NA, length(table$M)) # OUTPUT VECTOR
  
  for(i in 1:length(table$M)){
    S_curr <- RM_symm(table$M[i],f,ep)
    prop <- avgprop_real_components(evec_frame(S_curr))
    #print(prop)
    prop_vec[i] <- prop
  }
  cbind(table,prop_vec)
}

plot_M <- function(table, f){
  ggplot() + 
    geom_point(data = table, aes(x=M, y=prop_vec, color = prop_vec)) +
    labs(color = "EV Real", title = paste("f = ",f,sep="")) +
    scale_color_gradient(high="blue", low="red")
}

#### SYMM_EVC_PROPS

simulate_by_f_PROP2 <- function(f,M_max,ep_max,draws){
  M_vec <- sample(1:M_max, draws, replace = T)
  ep_vec <- sample(1:ep_max, draws, replace = F)
  table <- data.frame(M = M_vec, ep = rep(ep_vec,length(M_vec)))
  
  prop_vec <- rep(NA, length(table$M)) # OUTPUT VECTOR
  
  for(i in 1:length(table$M)){
    S_curr <- RM_symm(table$M[i],f,table$ep[i])
    prop <- avgprop_real_components(evec_frame(S_curr))
    #print(prop)
    prop_vec[i] <- prop
  }
  cbind(table,prop_vec)
}

props_M_ep <- function(table, f){
  ggplot() + 
    geom_point(data = table, aes(x=M, y=ep, color = prop_vec)) +
    labs(color = "EV Real", title = paste("f = ",f,sep="")) +
    scale_color_gradient(high="blue", low="red")
}

### SYMM_EVC

simulate_by_f_BOOL <- function(f,M_max,ep_max,draws){
  M_vec <- sample(1:M_max, draws, replace = T)
  ep_vec <- sample(1:ep_max, draws, replace = F)
  table <- data.frame(M = M_vec, ep = rep(ep_vec,length(M_vec)))
  
  bool_vec <- rep(NA, length(table$M))
  
  for(i in 1:length(table$M)){
    S_curr <- RM_symm(table$M[i],f,table$ep[i])
    bool_vec[i] <- check_real_eigenvectors(evec_frame(S_curr))
  }
  cbind(table,bool_vec)
}

plot_f_table <- function(table, f){
  ggplot() + 
    geom_point(data = table, aes(x=M, y=ep, color = factor(bool_vec))) +
    labs(color = "EV Real", title = paste("f = ",f,sep=""))
}

#=================================================================================#
#                           RANDOM SYMMETRIC MATRICES
#=================================================================================#

RM_symm <- function(M,f,ep){
  P <- matrix(data = runif(M**2, ep*(f-1), ep*f), nrow = M)
  .make_hermitian(P)
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

### FOUND IN LEGACY SYMM_EVC

unif_fpos <- function(M,f=T,ep){
  # unless specifically initialized, a random fraction will be chosen
  if(f){
    f <- runif(1,0,1)
    paste("f: ",f,sep="")
  }
  b <- f
  a <- (f-1)
  dist <- data.frame(x = runif(M**2, ep*a, ep*b))
  dist <- dist %>% mutate(x_neg = ifelse(x < 0,yes = 1, no = 0))
  dist
}