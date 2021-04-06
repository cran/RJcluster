RJ_mean_c = function(K, class, GG, GG_new)
{
  N = length(class)
  
  values = table(class)
  real_class = vector()
  for (i in 1:length(values))
  {
    real_class = c(real_class, rep(i, as.numeric(values[i])))
  }
  class = real_class

  # get the gamma values  
  gamma = get_gamma_c(class, K)
  gamma = gamma[unlist(lapply(gamma, length)) > 0]
  
  K = length(unique(class)) # redefine K 
  temp = get_mean_values_c(K, gamma, GG)
  mean_off = temp$mean_off
  mean_diag = temp$mean_diag
  
  # mean_diag     # real diagonals 
  # mean_off      # diagonals are homogeneous off diagonals
  if (K > 1)
  {
    mean_off = mean_off + t(mean_off) - diag(diag(mean_off))
  }
  
  # get mu values
  MU = get_mu_c(class, gamma, GG_new, N, K, mean_off, mean_diag)
  
  return(MU)
}


