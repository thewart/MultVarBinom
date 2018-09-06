# used_obs <- all_obs[all_obs[,length(Observation),by=FocalID][V1>=120,FocalID]]
# used_obs <- all_obs
library(standardize)
source("~/code/OrdRegMix/messyprep.R")

behaviors = c("SDB","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)", "Approach:initiate(partner)",
              "NonConAgg_give","NonConAgg_rec","contactAgg:direct'n(give)","contactAgg:direct'n(receive)")

used_obs <- all_obs[Group=="F" & SEX=="f"]
unique_obs <- used_obs[,.(FocalID,Year)] %>% unique()
std <- standardize(rep(1,nrow(unique_obs))~Year,unique_obs,family=binomial)
Xf <- model.matrix(lm(formula=std$formula,data=std$data))[,-1] %>% t()
Xr <- model.matrix( ~ 0 + FocalID,data=unique_obs) %>% t()
Y <- used_obs[,behaviors,with=F] %>% as.matrix() %>% t()
n <- used_obs[,length(Observation),by=.(FocalID,Year)]$V1
M <- length(n)

mvb <- stan_model("~/code/MultVarBinom/mvbinom.stan")
mvbn <- stan_model("~/code/MultVarBinom/mvbinom_mm_null.stan")
mvbv <- stan_model("~/code/MultVarBinom/mvbinom_mm_hv.stan")

standat <- list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y),X=Xf,P=nrow(Xf))

fit <- sampling(mvb,c(standat,f2_sigma0=1),
                iter=1000,chains=4,warmup=250,thin=3,include=F,cores=4)
fit0 <- sampling(mvb,c(standat,f2_sigma0=1e-16),
                iter=1000,chains=4,warmup=250,thin=3,include=F,cores=4)

standat <- list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y),X=Xf,P=nrow(Xf),Z=Xr,R=nrow(Xr))
fitn <- sampling(mvbn,standat,pars=c("f1_u"),
                 iter=1000,warmup=250,thin=3,chains=4,include=F,cores=4)
fitv <- sampling(mvbv,standat,pars=c("f1_u","f2_u","f2_sigma_raw"),
                 iter=1000,warmup=250,thin=3,chains=4,include=F,cores=4)
