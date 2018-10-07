make_yset <- function(D) return(do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t())

make_xset <- function(yset) return(apply(yset,2,interact))

calc_p <- function(f1,f2,yset=make_yset(nrow(f1)),xset=make_xset(yset)) { #f1 must be either D-length vector or D x N matrix
  eta <- t(yset) %*% f1 + t(xset) %*% f2
  return(apply(eta,2,function(x) exp(x)/sum(exp(x))))
}

calc_marg <- function(p,yset=make_yset(log2(nrow(p)))) {
  D <- log2(nrow(p))
  N <- ncol(p)

  ymu <- matrix(nrow=D,ncol=N)
  for (i in 1:D) ymu[i,] <- p[yset[i,]==1,] %>% colSums()
  
  return(ymu)
}

calc_marg_f <- function(f1,f2,yset=make_yset(nrow(f1))) {
  p <- calc_p(f1,f2,yset)
  return(calc_marg(p,yset))
}

calc_cor <- function(p,yset=make_yset(log2(nrow(p)))) {
  D <- log2(nrow(p))
  N <- ncol(p)
  k <- 1
  
  ycor <- matrix(nrow=D*(D-1)/2,ncol=N)
  for (i in 1:(D-1)) {
    for (j in (i+1):D) {
      p11 <- p[yset[i,]==1 & yset[j,]==1,] 
      p00 <- p[yset[i,]==-1 & yset[j,]==-1,] 
      p01 <- p[yset[i,]==-1 & yset[j,]==1,] 
      p10 <- p[yset[i,]==1 & yset[j,]==-1,] 
      if (D>2) {
        if (N>1) {
          p11 <- colSums(p11)
          p00 <- colSums(p00)
          p01 <- colSums(p01)
          p10 <- colSums(p10)
        } else {
          p11 <- sum(p11)
          p00 <- sum(p00)
          p01 <- sum(p01)
          p10 <- sum(p10)
        }
      }
      ycor[k,] <- (p11*p00 - p10*p01)/sqrt((p11+p10)*(p00+p01)*(p11+p01)*(p00+p10))
      k <- k + 1
    }
  }
  return(ycor)
}

calc_cor_f <- function(f1,f2,yset=make_yset(nrow(f1))) {
  p <- calc_p(f1,f2,yset)
  return(calc_cor(p,yset))
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

calc_par_f <- function(FUN,f1,f2,cores=4) {
  n <- dim(f1)[1]
  m <- dim(f1)[3]
  
  f1l <- f2l <- list()
  for (i in 1:n) f1l[[i]] <- f1[i,,]
  for (i in 1:n) f2l[[i]] <- f2[i,,]
  
  out <- mcmapply(FUN,f1l,f2l,mc.cores=cores)
  
  d <- nrow(out)/m
  return(array(out,dim=c(d,m,n)))
}

calc_par <- function(FUN,p,cores=4) { #p dim: c(states, observations, samples)
  n <- dim(p)[3] #number of samples
  m <- dim(p)[2] #number of animals
  
  pl <- list()
  for (i in 1:n) pl[[i]] <- p[,,i]
  
  out <- mclapply(pl,FUN,mc.cores=cores)
  
}

make_pltdat <- function(dat,label=NA) {
  return(data.table(label=label,std=apply(dat,1,median),
                    ubi=apply(dat,1,quantile,0.1),ubo=apply(dat,1,quantile,0.025),
                    lbi=apply(dat,1,quantile,0.9),lbo=apply(dat,1,quantile,0.975)))
}

crossent <- function(p,q) {
  s <- p*log(q)
  return( -colSums(p*log(q)) )
}

match_marg <- function(p,yset=make_yset(log2(nrow(p)))) {
  yset[yset==-1] <- 0
  mup <- calc_marg(p,yset)
  return( exp(t(yset) %*% log(mup) + t(!yset) %*% log(1-mup)) )
}

calc_repeat <- function(p,mu=rowMeans(p)) {
  D0 <- crossent(p,mu) %>% mean()
  D <- crossent(p,p) %>% mean()
  return( (D0-D)/D0 )
}

calc_repeat_marg <- function(p) {
  D0 <- crossent(p,rowMeans(p)) %>% mean()
  D <- crossent(p,match_marg(p)) %>% mean()
  return( (D0-D)/D0 )
}
