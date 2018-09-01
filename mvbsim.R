library(rstan)

source("~/code/OrdRegMix/messyprep.R")

behaviors = c("SDB","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)", "Approach:initiate(partner)",
              "NonConAgg_give","NonConAgg_rec","contactAgg:direct'n(give)","contactAgg:direct'n(receive)")
D <- length(behaviors)
used_obs <- all_obs[Group=="F" & SEX=="f"]
Y <- used_obs[,behaviors,with=F] %>% as.matrix() %>% t()
N <- ncol(Y)
n <- used_obs[,length(Observation),by=FocalID]$V1
M <- length(n)

yset <- do.call(expand.grid,rep(list(c(-1,1)),D)) %>% t()
xset <- apply(yset,2,interact)

mvbn <- stan_model("~/code/MultVarBinom/mvbinom_mm_null.stan")
initfit <- sampling(mvbn,list(Y=Y,N=N,D=D,D2=2^D,M=M,n=n),
                chain=4,cores=4,iter=650,thin=2,warmup=150,pars=c("f1_raw","f2_raw","f2_sigma_raw"),include=F)

tau_f1_sigma <- 0.5
tau_f2_sigma <- 0.5

f1 <- init$f1 + f1_sigma*matrix(rnorm(D*M),nrow=D)
f2 <- init$f2 + f2_sigma*matrix(rnorm(M*D*(D-1)/2),nrow=D*(D-1)/2)

p <- t(yset) %*% f1 + t(xset) %*% f2
p <- apply(p,2,function(x) exp(x)/sum(exp(x)))

ycat <- matrix(nrow=2^D,ncol=0)
for (i in 1:M) {
  ycat <- cbind(ycat,rmultinom(n[i],1,p[,i]))
}
y <- yset[,apply(ycat==1,2,which)]
y[y==-1] <- 0

mvb <- stan_model("~/code/MultVarBinom/mvbinom.stan")
init <- optimizing(mvb,list(Y=y,N=N,D=D,D2=2^D),as_vector=F)$par

mvbm <- stan_model("~/code/MultVarBinom/mvbinom_mm.stan")
fit <- sampling(mvbm,list(Y=y,N=N,D=D,D2=2^D,M=M,n=n),
                chain=3,cores=3,iter=1000,init=rep(list(init),3),pars=c("f1_raw","f2_raw"),include=F)
mvbv <- stan_model("~/code/MultVarBinom/mvbinom_mm_hv_c.stan")
fit <- sampling(mvbv,list(Y=y,N=N,D=D,D2=2^D,M=M,n=n),
                chain=3,cores=3,iter=400,warmup=150,pars=c("f1_raw","f2_raw","f2_sigma_raw"),include=F)


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
