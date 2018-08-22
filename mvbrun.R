# used_obs <- all_obs[all_obs[,length(Observation),by=FocalID][V1>=120,FocalID]]
# used_obs <- all_obs
used_obs <- all_obs[Group=="F"]
Y <- used_obs[,.(GroomGET,GroomGIVE,`Approach:initiate(focal)`,passcont)] %>% as.matrix() %>% t()
n <- used_obs[,length(Observation),by=FocalID]$V1
M <- length(n)

mvb <- stan_model("~/code/multvarbinom/mvbinom.stan")
fit <- sampling(mvb,list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n),iter=500,chains=3,
                pars=c("f1_mu","f1_sigma","f2_mu","f2_sigma","f1","f2"),cores=4)

