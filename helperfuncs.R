calc_p <- function(f1,f2) { #f1 must be either D-length vector or D x N matrix
  D <- length(f1)
  yset <- do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t()
  xset <- apply(yset,2,interact)
  
  eta <- t(yset) %*% f1 + t(xset) %*% f2
  return(apply(eta,2,function(x) exp(x)/sum(exp(x))))
}

calc_marg_cov <- function(f1,f2) {
  D <- length(f1)
  yset <- do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t()
}
