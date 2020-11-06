
# This script includes functions that are metrics for analyzing eigenvectors.

# The inputs for these metric functions have been standardized to the eigen_frame object returned 
# by the eigen_frame() function in eigenvectors.R.

#=================================================================================#
#                               EIGENVECTOR METRICS
#=================================================================================#

L_eigenplot <- function(P){
  eigen_plot(t(P), mat_type = "Transpose")
}

R_eigenplot <- function(P){
  eigen_plot(P, mat_type = "Original")
  }

check_real_eigenvectors <- function(eigen_frame){
  prop <- prop_real_rows(eigen_frame)
  real_prop <- round(sum(prop$is_real)/length(prop$is_real),10)
  ifelse(real_prop == 1, T, F)
}

prop_real_rows <- function(eigen_frame){
  M <- sqrt(length(eigen_frame$Re)) #obtain number of rows
  eigen_frame %>%
    mutate(Im_0 = case_when(Im == 0 ~ T, Im != 0 ~ F)) %>%
    group_by(row_i) %>%
    summarize(prop_reals = round(sum(Im_0)/M, 4)) %>%
    mutate(is_real = case_when(prop_reals == 1 ~ T, prop_reals < 1 ~ F))
}

norms <- function(eigen_frame){
  #returns maximum norm by component, maximum norm by row,
  #average norm of component, row
  1
}

#=================================================================================#
#                         AGGREGATE/SUMMARY EIGENMETRICS
#=================================================================================#

rmt_summary <- function(P){
  evs <- eigen(P)[1]$values
  RMThreshold::rm.spacing.distribution(evs)
  RMThreshold::rm.ev.density(evs)
}

eigen_summary <- function(eigen_frame, loud = F){
  prop <- prop_real_rows(eigen_frame)
  if(loud){print(prop)}
  real_prop <- round(sum(prop$is_real)/length(prop$is_real),10)
  paste("Proportion of real-valued rows: ",real_prop,sep="")
}
