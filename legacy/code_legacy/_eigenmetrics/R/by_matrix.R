
# This script includes functions that are matrix specific functions for analyzing eigenvectors of a given type of matrix.

#=================================================================================#
#                                     MATRIXES
#=================================================================================#

symm_stoch <- function(M, bool_plot=T, bool_loud=F){
  P <- rand_M_symm_stoch(M, row_fn = r_zeros)
  eigen_summary(eigen_frame(P))
  rmt_summary(P)
  if(bool_plot){eigen_plot(P, loud = bool_loud, "Symmetric Stochastic")}
}

symm_norm <- function(M, bool_plot=T, bool_loud=F){
  P <- rand_M_symm_norm(M, mu = 0, sd = 1)
  if(bool_plot){eigen_plot(P, loud = bool_loud, "Normal Symmetric")}
  eigen_summary(eigen_frame(P))
  rmt_summary(P)
}

tridiag <- function(M, bool_plot=T, bool_loud=F){
  P <- rand_M_trid(M)
  if(bool_plot){eigen_plot(P, loud = bool_loud, "Tridiagonal")}
  eigen_summary(eigen_frame(P))
  rmt_summary(P)
}