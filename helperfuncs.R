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

interact_string <- function(y) {
  D = length(y)
  x <- vector(length=D*(D-1)/2)
  k <- 1
  
  for (i in 1:(D-1)) {
    for (j in (i+1):D) {
      x[k] <- paste0(y[i],":",y[j])
      k <- k + 1;
    }
  }
  return(x);
}

calc_par_f <- function(FUN,f1,f2) {
  n <- dim(f1)[1]
  m <- dim(f1)[3]
  
  f1l <- f2l <- list()
  for (i in 1:n) f1l[[i]] <- f1[i,,]
  for (i in 1:n) f2l[[i]] <- f2[i,,]
  
  out <- mapply(FUN,f1l,f2l)
  
  d <- nrow(out)/m
  return(array(out,dim=c(d,m,n)))
}

calc_par <- function(FUN,p,return_array=T) { #p dim: c(states, observations, samples)
  n <- dim(p)[3] #number of samples
  m <- dim(p)[2] #number of animals
  
  pl <- list()
  for (i in 1:n) pl[[i]] <- p[,,i]
  
  out <- lapply(pl,FUN)
  if (return_array) out <- unlist(out) %>% array(dim = c(dim(out[[1]]),length(out)))
}

calc_par_f <- function(FUN,f1,f2) {
  n <- dim(f1)[1]
  m <- dim(f1)[3]
  f1l <- f2l <- list()
  for (i in 1:n) f1l[[i]] <- f1[i,,]
  for (i in 1:n) f2l[[i]] <- f2[i,,]
  out <- mapply(FUN,f1l,f2l)
  d <- nrow(out)/m
  return(array(out,dim=c(d,m,n)))
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

calc_mi <- function(p,normalize=T) {
  mu <- rowMeans(p)
  D0 <- -sum(mu*log(mu))
  D <- crossent(p,p) %>% mean()
  mi <- D0-D
  if (normalize) mi <- mi/D0
  return(mi)
}

calc_mi_marg <- function(p,normalize=T) {
  pmarg <- match_marg(p)
  D0 <- -sum(rowMeans(p)*log(rowMeans(pmarg)))
  D <- crossent(p,match_marg(p)) %>% mean()
  mi <- D0-D
  if (normalize) mi <- mi/D0
  return(D0-D)
}

calc_popkldiv <- function(p,rank=F,forplt=T) {
  klmat <- apply(p,3,function(x) crossent(rowMeans(x),x)+sum(rowMeans(x)*log(rowMeans(x))))
  if (rank) klmat <- apply(klmat,2,frankv,order=-1)/nrow(klmat)
  if (forplt) klmat <- mat2dat(klmat)
  return(klmat)
}

calc_p_mm <- function(fit) {
  samps <- extract(fit,pars=c("f1_mu","f1_u","f2_mu","f2_u"))
  n <- dim(samps$f1_u)[1]
  d <- dim(samps$f1_u)[2]
  m <- dim(samps$f1_u)[3]
  f1 <- array(dim=dim(samps$f1_u))
  f2 <- array(dim=dim(samps$f2_u))
  for (i in 1:n) {
    f1[i,,] <- samps$f1_mu[i,] + samps$f1_u[i,,]
    f2[i,,] <- samps$f2_mu[i,] + samps$f2_u[i,,]
  }
  
  return(calc_par_f(calc_p,f1,f2))
}

calc_p_f1 <- function(fit) {
  samps <- extract(fit,pars=c("f1_mu","f1_u","f2"))
  n <- dim(samps$f1_u)[1]
  d <- dim(samps$f1_u)[2]
  m <- dim(samps$f1_u)[3]
  f1 <- array(dim=dim(samps$f1_u))
  for (i in 1:n) f1[i,,] <- samps$f1_mu[i,] + samps$f1_u[i,,]
  f2 <- array(samps$f2,dim=c(dim(samps$f2),m))
  
  return(calc_par_f(calc_p,f1,f2))
}

mat2dat <- function(X,dim=1) {
  return(data.table(kldiv=apply(X,dim,mean),lb=apply(X,dim,quantile,probs=0.1),llb=apply(X,dim,quantile,probs=0.025),
                    ub=apply(X,dim,quantile,probs=0.9),uub=apply(X,dim,quantile,probs=0.975)))
}