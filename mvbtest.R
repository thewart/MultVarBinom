library(rstan)
N <- 1000
M <- 10
n <- rep(N/M,M)
p <- rnorm(4)
p <- exp(p)/sum(exp(p))

yset <- t(expand.grid(c(0,1),c(0,1)))
ycat <- rmultinom(N,1,p)
y <- yset[,apply(ycat==1,2,which)]


mvb <- stan_model("~/code/multvarbinom/mvbinom.stan")
fit <- sampling(mvb,list(Y=y,N=N,D=nrow(y),M=M,n=n),iter=100)
fit <- optimizing(mvb,list(Y=y,N=N,D=nrow(y),M=M,n=n),as.vector=F)

i <- 2
phat <- c(0,fit$par$f1[i,1],fit$par$f1[i,2],sum(fit$par$f1[i,])+fit$par$f2[i])
phat <- exp(phat)/sum(exp(phat))
phat - rowMeans(ycat[,((i-1)*100+1):(i*100)])


