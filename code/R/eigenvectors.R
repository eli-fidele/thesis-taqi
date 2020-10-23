
# This script includes functions that help extract and plot eigenvalues of transition matrices.

#=================================================================================#
#                               EIGENVECTOR FUNCTIONS
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

# Plots the eigenvalues of a given matrix P
eigen_plot <- function(P, loud=F, mat_type=""){
  M <- length(P[1,])
  # Create eigenvalue dataframe using eigen_frame(P) function defined above
  eigen_frame_P <- eigen_frame(P)
  # Print eigenvalues and corrosponding matrix
  if(loud == T){ 
    if(M <= 10){print(eigen_frame_P)}            
    if(M <= 6){print(P)}
  }
  # Plot parameters
  r <- 1
  r_ep <- r + 0.5
  # Plot
  ggplot(eigen_frame_P) + 
    geom_point(aes(x = Re, y = Im, color = factor(row_i))) + 
    labs(x = "Re", y = "Im", title = paste("Eigenvectors: ",mat_type," Matrix",sep = "")) +
    xlim(-r_ep,r_ep) + ylim(-r,r) +
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r), color = "steelblue") +
    scale_color_discrete(name = "Row_i") +
    scale_x_continuous(name = "Re(r_i)") +
    scale_y_continuous(name = "Im(r_i)") +
    coord_fixed(ratio = 1)
}