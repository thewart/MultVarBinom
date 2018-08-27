library(rstan)
N <- 10000
# M <- 10
# n <- rep(N/M,M)
D <- 10
yset <- do.call(expand.grid,rep(list(c(0,1)),D)) %>% t()
xset <- apply(yset,2,interact)
f1 <- rnorm(D)*0.5
f2 <- rnorm(D*(D-1)/2)*0.1

p <- (t(f1) %*% yset + t(f2) %*% xset) %>% as.vector()
Z <- sum(exp(p))
p <- exp(p)/Z

ycat <- rmultinom(N,1,p)
y <- yset[,apply(ycat==1,2,which)]

mvb <- stan_model("~/code/MultVarBinom/mvbinom.stan")
fit <- sampling(mvb,list(Y=y,N=N,D=D),chain=1,iter=100)
fit <- optimizing(mvb,list(Y=y,N=N,D=D),as_vector=F,verbose=T)

mvbs <- stan_model("~/code/MultVarBinom/mvbinom_store.stan")
fit <- sampling(mvbs,list(Y=y,N=N,D=D,D2=2^D),chain=1,iter=100)


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
