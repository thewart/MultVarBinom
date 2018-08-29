# used_obs <- all_obs[all_obs[,length(Observation),by=FocalID][V1>=120,FocalID]]
# used_obs <- all_obs
used_obs <- all_obs[Group=="F" & SEX=="f"]
Y <- used_obs[,behaviors,with=F] %>% as.matrix() %>% t()
n <- used_obs[,length(Observation),by=FocalID]$V1
M <- length(n)

mvb <- stan_model("~/code/MultVarBinom/mvbinom_mm.stan")
fitm <- sampling(mvb,list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y)),iter=1000,chains=5,
                pars=c("f1_raw","f2_raw"),include=F,cores=5)

mvbv <- stan_model("~/code/MultVarBinom/mvbinom_mm_hv.stan")
fitv <- sampling(mvbv,list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y)),iter=1000,chains=5,
                pars=c("f1_raw","f2_raw","f2_sigma_raw"),include=F,cores=5)
