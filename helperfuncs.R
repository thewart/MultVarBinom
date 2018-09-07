make_yset <- function(D) return(do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t())

make_xset <- function(yset) return(apply(yset,2,interact))

calc_p <- function(f1,f2,yset=make_yset(nrow(f1)),xset=make_xset(yset)) { #f1 must be either D-length vector or D x N matrix
  eta <- t(yset) %*% f1 + t(xset) %*% f2
  return(apply(eta,2,function(x) exp(x)/sum(exp(x))))
}

calc_marg <- function(f1,f2,yset=make_yset(nrow(f1)),p=calc_p(f1,f2,yset)) {
  D <- nrow(f1)
  N <- ncol(f1)

  ymu <- matrix(nrow=D,ncol=N)
  for (i in 1:D) ymu[i,] <- p[yset[i,]==1,] %>% colSums()
  
  return(ymu)
}

calc_cor <- function(f1,f2,yset=make_yset(nrow(f1)),p=calc_p(f1,f2,yset)) {
  D <- nrow(f1)
  N <- ncol(f1)
  k <- 1
  
  ycor <- matrix(nrow=D*(D-1)/2,ncol=N)
  for (i in 1:(D-1)) {
    for (j in (i+1):D) {
      p11 <- p[yset[i,]==1 & yset[j,]==1,] %>% colSums()
      p00 <- p[yset[i,]==-1 & yset[j,]==-1,] %>% colSums()
      if (i != j) {
        p01 <- p[yset[i,]==-1 & yset[j,]==1,] %>% colSums()
        p10 <- p[yset[i,]==1 & yset[j,]==-1,] %>% colSums()
      } else {
        p01 <- p10 <- rep(0,N)
      }
      ycor[k,] <- (p11*p00 - p10*p01)/sqrt((p11+p10)*(p00+p01)*(p11+p01)*(p00+p10))
      k <- k + 1
    }
  }
  return(ycor)
}

interact <- function(y) {
  D = length(y)
  x <- vector(length=D*(D-1)/2)
  k <- 1
  
  for (i in 1:(D-1)) {
    for (j in (i+1):D) {
      x[k] <- y[i]*y[j]
      k <- k + 1;
    }
  }
  return(x);
}

