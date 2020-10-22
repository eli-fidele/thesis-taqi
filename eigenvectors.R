
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
eigen_plot <- function(P, loud){
  # Create eigenvalue dataframe using eigen_frame(P) function defined above
  eigen_frame_P <- eigen_frame(P)
  # Print eigenvalues and corrosponding matrix
  if(loud == T){ # Check if it possible to set loud == F as a default
    print(eigen_frame_P)            
    print(P)
  }
  # Plot
  ggplot(eigen_frame_P) + 
    geom_point(aes(x = Re, y = Im, color = factor(row_i))) + 
    labs(x = "Re", y = "Im", title = "Distribution of Eigenvectors in C") +
    xlim(-(r+ep),r+ep) + ylim(-r,r) +
    ggforce::geom_circle(aes(x0 = 0, y0 = 0, r = r), color = "blue") +
    scale_color_discrete(name = "Row_i") +
    scale_x_continuous(name = "Re(r_i)") +
    scale_y_continuous(name = "Im(r_i)") +
    coord_fixed(ratio = 1)
}