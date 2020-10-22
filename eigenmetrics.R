
# This script includes functions that are metrics for analyzing eigenvectors.

# The inputs for these metric functions have been standardized to the eigen_frame object returned 
# by the eigen_frame() function in eigenvectors.R.

#=================================================================================#
#                               EIGENVECTOR METRICS
#=================================================================================#

eigen_summary <- function(eigen_frame){
  prop <- prop_real_rows(eigen_frame)
  print(prop)
  real_prop <- round(sum(prop$is_real)/length(prop$is_real),10)
  paste("Proportion of real-valued rows: ",real_prop,sep="")
}

prop_real_rows <- function(eigen_frame){
  M <- sqrt(length(eigen_frame$Re)) #obtain number of rows
  eigen_frame %>%
    mutate(Im_0 = case_when(Im == 0 ~ T, Im != 0 ~ F)) %>%
    group_by(row_i) %>%
    summarize(prop_reals = round(sum(Im_0)/M, 4)) %>%
    mutate(is_real = case_when(prop_reals == 1 ~ T, prop_reals < 1 ~ F))
}