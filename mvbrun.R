# used_obs <- all_obs[all_obs[,length(Observation),by=FocalID][V1>=120,FocalID]]
# used_obs <- all_obs
library(standardize)
library(lme)
source("~/code/OrdRegMix/messyprep.R")

behaviors = c("SDB","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)", "Approach:initiate(partner)",
              "NonConAgg_give","NonConAgg_rec","contactAgg:direct'n(give)","contactAgg:direct'n(receive)")

used_obs <- all_obs[Group=="F" & SEX=="f"]
std <- standardize(SDB~Year,used_obs,family=binomial)
Xf <- model.matrix(lm(formula=std$formula,data=std$data))[,-1]
Y <- used_obs[,behaviors,with=F] %>% as.matrix() %>% t()
n <- used_obs[,length(Observation),by=FocalID]$V1
M <- length(n)

mvb <- stan_model("~/code/MultVarBinom/mvbinom.stan")
mvbn <- stan_model("~/code/MultVarBinom/mvbinom_mm_null.stan")
mvbv <- stan_model("~/code/MultVarBinom/mvbinom_mm_hv.stan")

fit <- sampling(mvb,list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y)),
                 iter=1000,chains=4,warmup=250,thin=3,include=F,cores=4)

fitn <- sampling(mvbn,list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y)),pars=c("f1_raw","f2_raw"),
                 iter=1000,warmup=250,thin=3,chains=4,include=F,cores=4)

fitv <- sampling(mvbv,list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y),X=t(Xf),P=ncol(Xf)),pars=c("f1_raw","f2_raw","f2_sigma_raw"),
                 iter=1000,warmup=250,thin=3,chains=4,include=F,cores=4)
