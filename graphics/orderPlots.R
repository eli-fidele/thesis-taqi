# ===================================================
#                       Setup
# ===================================================
# Real histogram
plot_h_Re <- function(ens_spec, N, ordering){
  ens_spec %>% 
    mutate(Order = as.factor(Order)) %>% 
    ggplot(mapping = aes(x = Re, fill = Order)) +
    geom_histogram(alpha = alpha, bins = bins, inherit.aes = T) +
    scale_fill_manual(values = wes_palette("Zissou1", N, type = "continuous"))+
    labs(fill = "Order", 
         title = paste(ordering,"-Ordered Spectrum of a ", 
                       makeTitle(N, "Symmetric Normal Matrix Ensemble"), sep = ""))
}
# Real density
plot_d_Re <- function(ens_spec, N, ordering){
  ens_spec <- ens_spec %>% 
    mutate(Order = as.factor(Order)) %>%
    ggplot(mapping = aes(x = Re, fill = Order)) +
    geom_density(alpha = alpha, bins = bins, inherit.aes = T)  +
    labs(fill = "Order", 
         title = paste(ordering,"-Ordered Spectrum of a ", 
                       makeTitle(N, "Symmetric Normal Matrix Ensemble"), sep = ""))
}
# Norm density
plot_d_Norm <- function(ens_spec, N, ordering){
  ens_spec %>%
    mutate(Order = as.factor(Order)) %>%
    ggplot(mapping = aes(x = Norm, group = Order, fill = Order)) +
    geom_density(alpha = alpha, kernel = "gaussian") +
    labs(fill = "Order", 
         title = paste(ordering,"-Ordered Spectrum of a ", 
                       makeTitle(N, "Symmetric Normal Matrix Ensemble"), sep = ""))
}
# Norm density
plot_d_Norm1 <- function(ens_spec, N){
  ens_spec %>% 
    mutate(Order = as.factor(Order)) %>%
    ggplot(mapping = aes(x = Norm, fill = Order)) +
    geom_density(alpha = alpha, bins = bins, inherit.aes = T) +
    scale_fill_manual(values = wes_palette("Zissou1", N, type = "continuous"))+
    labs(fill = "Order", 
         title = paste("Norm-ordered Spectrum of a", 
                       makeTitle(N, "Symmetric Normal Matrix Ensemble")))
}
plot_d_Norm2 <- function(ens_spec, N){
  ens_spec %>% 
    mutate(Order = as.factor(Order)) %>% 
    ggplot(mapping = aes(x = Norm, fill = Order)) +
    geom_density(alpha = alpha, bins = bins, inherit.aes = T) +
    scale_fill_manual(values = wes_palette("Zissou1", N, type = "continuous"))
}