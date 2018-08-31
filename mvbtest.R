library(rstan)
N <- 10000
M <- 50
n <- rep(N/M,M)
D <- 10
yset <- do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t()
xset <- apply(yset,2,interact)
f1 <- rnorm(D)*0.5
f2 <- rnorm(D*(D-1)/2)*0.25

p <- (t(f1) %*% yset + t(f2) %*% xset) %>% as.vector()
Z <- sum(exp(p))
p <- exp(p)/Z

ycat <- rmultinom(N,1,p)
y <- yset[,apply(ycat==1,2,which)]
y[y==-1] <- 0

mvb <- stan_model("~/code/MultVarBinom/mvbinom.stan")
init <- optimizing(mvb,list(Y=y,N=N,D=D,D2=2^D),as_vector=F)$par

mvbm <- stan_model("~/code/MultVarBinom/mvbinom_mm.stan")
fit <- sampling(mvbm,list(Y=y,N=N,D=D,D2=2^D,M=M,n=n),
                chain=3,cores=3,iter=1000,init=rep(list(init),3),pars=c("f1_raw","f2_raw"),include=F)
mvbv <- stan_model("~/code/MultVarBinom/mvbinom_mm_hv.stan")
fit <- sampling(mvbv,list(Y=y,N=N,D=D,D2=2^D,M=M,n=n),
                chain=3,cores=3,iter=1000,init=rep(list(init),3),pars=c("f1_raw","f2_raw"),include=F)


interact <- function(y) {
  D <- length(y);
  x <- vector(length=D*(D-1)/2)
  k <- 1
  for (i in 1:(D-1)) {
    for (j in (i+1):D) {
      x[k] = y[i]*y[j];
      k = k + 1;
    }
  }
  return(x)
}
