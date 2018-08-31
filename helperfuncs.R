calc_p <- function(f1,f2) { #f1 must be either D-length vector or D x N matrix
  D <- length(f1)
  yset <- do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t()
  xset <- apply(yset,2,interact)
  
  eta <- t(yset) %*% f1 + t(xset) %*% f2
  return(apply(eta,2,function(x) exp(x)/sum(exp(x))))
}

calc_marg_cov <- function(f1,f2) {
  D <- nrow(f1)
  N <- ncol(f1)
  yset <- do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t()
  
  p <- calc_p(f1,f2)
  
  ymu <- matrix(nrow=D,ncol=N)
  for (i in 1:D) ymu[d,] <- p[yset[i,]==1,] %>% rowSums()
  
  ycor <- array(dim = c(D,D,N))
  for (i in 1:(D-1)) {
    for (j in i:D) {
      p11 <- p[yset[i,]==1 & yset[j,]==1,] %>% rowSums()
      p10 <- p[yset[i,]==1 & yset[j,]==-1,] %>% rowSums()
      p01 <- p[yset[i,]==-1 & yset[j,]==1,] %>% rowSums()
      p00 <- p[yset[i,]==-1 & yset[j,]==-1,] %>% rowSums()
      ycor[j,i,] <- ycor[i,j,] <- (p11*p00 - p10*p01)/sqrt((p11+p10)*(p00+p01)*(p11+p01)*(p00+p10))
    }
  }
}
