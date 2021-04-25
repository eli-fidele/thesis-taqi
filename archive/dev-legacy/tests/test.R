
# Test proportion of numerically-valued (valid) symmetric stochastic matrices
stoch.VALID <- function(iter, N, symm = F){
  mean(replicate(n = iter, .pos_entries(RM_stoch(N, symm), zero = T) == 1, simplify = T))
}

# Test row-stochasticity of symmetric stochastic matrices
stoch.STOCH <- function(iter, N, symm = F){
  mean(replicate(n = 100, .isStochastic(RM_stoch(3, symm = T)), simplify = T))
}

norm.SYMM <- function(iter, N, cplx = F){
  mean(map_lgl(RME_norm(N, symm = T, cplx = cplx, herm = cplx, size = iter), .isHermitian))
}

norm.ENS1 <- function(iter, N, cplx = F){
  ens <- RME_norm(N, symm = T, cplx = cplx, herm = cplx, size = iter)
  ens
}

ens.SYMM <- function(ens){
  mean(map_lgl(ens, .isHermitian))
}