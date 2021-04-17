
# Plot settings
device = "pdf"
width = 20
height = 10
units = "cm"

library(RMAT)
set.seed(23)
# Generate random matrix and ensemble
P <- RM_norm(N = 5)
ens <- RME_norm(N = 5, size = 10)
# Get spectra
P_spec <- P %>% spectrum()
ens_spec <- ens %>% spectrum()
# Spectum plots
P_plot <- P_spec %>% spectrum.scatterplot()
ens_plot <- ens_spec %>% spectrum.scatterplot()
# Comparison
P_plot + ens_plot
ggsave("2-1-2_comparison.pdf")
