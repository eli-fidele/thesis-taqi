# ===================================================
#                       Setup
# ===================================================

# Labels a generic array with a value column
.label_array <- function(array, value, label, size){
  label_column <- data.frame(x = rep(value, size))
  colnames(label_column) <- label
  cbind(array, label_column)
}
# Note: Add this feature to spectrum.scatterplot()!
# Creates a matrix title given the dimension and classification of the matrix (ex. 20x20 Unif(0,1))
.make_matrix_title <- function(N, class){
  dim_str <- paste(N,"x",N, sep = "")
  title_str <- paste(dim_str, class, sep = " ")
  return(title_str)
}
# Create a shorter alias of this function
makeTitle <- .make_matrix_title 
# ===================================================
#                       Plots
# ===================================================
# Real histogram
plot_h_Re <- function(ens_spec, ordering, class = "Symmetric Normal Matrix Ensemble"){
  # Get matrix dimension
  N <- which.max(ens_spec$Order)
  # Return spectrum plot
  ens_spec %>% 
    mutate(Order = as.factor(Order)) %>% 
    ggplot(mapping = aes(x = Re, fill = Order)) +
    geom_histogram(alpha = alpha, bins = bins, inherit.aes = T) +
    scale_fill_manual(values = wes_palette("Zissou1", N, type = "continuous"))+
    labs(fill = "Order", 
         title = paste(ordering,"-Ordered Spectrum of a ", 
                       makeTitle(N, class), sep = ""))
}
# Real density
plot_d_Re <- function(ens_spec, ordering, class = "Symmetric Normal Matrix Ensemble"){
  # Get matrix dimension
  N <- which.max(ens_spec$Order)
  # Return spectrum plot
  ens_spec %>% 
    mutate(Order = as.factor(Order)) %>%
    ggplot(mapping = aes(x = Re, fill = Order)) +
    geom_density(alpha = alpha, bins = bins, inherit.aes = T)  +
    labs(fill = "Order", 
         title = paste(ordering,"-Ordered Spectrum of a ", 
                       makeTitle(N, class), sep = ""))
}
# Norm density
plot_d_Norm <- function(ens_spec, ordering, class = "Symmetric Normal Matrix Ensemble"){
  # Get matrix dimension
  N <- which.max(ens_spec$Order)
  # Return spectrum plot
  ens_spec %>%
    mutate(Order = as.factor(Order)) %>%
    ggplot(mapping = aes(x = Norm, group = Order, fill = Order)) +
    geom_density(alpha = alpha, kernel = "gaussian") +
    labs(fill = "Order", 
         title = paste(ordering,"-Ordered Spectrum of a ", makeTitle(N, class), sep = ""))
}
# Norm density
# plot_d_Norm1 <- function(ens_spec, N){
#   ens_spec %>% 
#     mutate(Order = as.factor(Order)) %>%
#     ggplot(mapping = aes(x = Norm, fill = Order)) +
#     geom_density(alpha = alpha, bins = bins, inherit.aes = T) +
#     scale_fill_manual(values = wes_palette("Zissou1", N, type = "continuous"))+
#     labs(fill = "Order", 
#          title = paste("Norm-ordered Spectrum of a", 
#                        makeTitle(N, "Symmetric Normal Matrix Ensemble")))
# }
# plot_d_Norm2 <- function(ens_spec, N){
#   ens_spec %>% 
#     mutate(Order = as.factor(Order)) %>% 
#     ggplot(mapping = aes(x = Norm, fill = Order)) +
#     geom_density(alpha = alpha, bins = bins, inherit.aes = T) +
#     scale_fill_manual(values = wes_palette("Zissou1", N, type = "continuous"))
# }