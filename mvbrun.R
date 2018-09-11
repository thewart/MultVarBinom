# used_obs <- all_obs[all_obs[,length(Observation),by=FocalID][V1>=120,FocalID]]
# used_obs <- all_obs
library(standardize)
library(shinystan) ##somehow prevents crashing on parallel chains in R 3.5?
source("~/code/OrdRegMix/messyprep.R")

behaviors = c("SDB","GroomGIVE", "GroomGET","passcont","Approach:initiate(focal)", "Approach:initiate(partner)",
              "NonConAgg_give","NonConAgg_rec","contactAgg:direct'n(give)","contactAgg:direct'n(receive)")

used_obs <- all_obs[Group=="F" & SEX=="f"]
unique_obs <- used_obs[,.(FocalID,Year)] %>% unique()
std <- standardize(rep(1,nrow(unique_obs))~Year,unique_obs,family=binomial)
Xf <- model.matrix(lm(formula=std$formula,data=std$data))[,-1] %>% t()
Xr <- model.matrix( ~ 0 + FocalID,data=unique_obs) %>% as.matrix() %>% t()
Y <- used_obs[,behaviors,with=F] %>% as.matrix() %>% t()
n <- used_obs[,length(Observation),by=.(FocalID,Year)]$V1
M <- length(n)

mvb <- stan_model("~/code/MultVarBinom/mvbinom.stan")
mvb_mm <- stan_model("~/code/MultVarBinom/mvbinom_mm.stan")
mvb_mm_f1 <- stan_model("~/code/MultVarBinom/mvbinom_mm_f1.stan")
mvb_mm_f2 <- stan_model("~/code/MultVarBinom/mvbinom_mm_f2.stan")

standat <- list(Y=Y,N=ncol(Y),D=nrow(Y),M=length(n),n=n,D2=2^nrow(Y),X=Xf,P=nrow(Xf))

fit0 <- sampling(mvb,standat,
                 iter=700,warmup=200,thin=2,chains=4,cores=4)

standat <- c(standat,list(Z=Xr,R=nrow(Xr)))
fit_mm_f2 <- sampling(mvb_mm_f2,standat,pars=c("f2_u_raw","f2_sigma_raw"),
                      iter=700,warmup=200,thin=2,chains=4,include=F,cores=4)
fit_mm_f1 <- sampling(mvb_mm_f1,standat,pars=c("f1_u_raw"),
                      iter=700,warmup=200,thin=2,chains=4,include=F,cores=4)
fit_mm <- sampling(mvb_mm,standat,pars=c("f1_u_raw","f2_u_raw","f2_sigma_raw"),
                   iter=700,warmup=200,thin=2,chains=4,include=F,cores=4)
