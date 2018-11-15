# used_obs <- all_obs[all_obs[,length(Observation),by=FocalID][V1>=120,FocalID]]
# used_obs <- all_obs
library(standardize)
library(shinystan) ##somehow prevents crashing on parallel chains in R 3.5?
modonsex <- T
modonkin <- T
modonrank <- F
source("~/code/OrdRegMix/collectcleanfocaldata.R")
source("~/code/OrdRegMix/loadcovariates.R")
source("~/code/OrdRegMix/tabulatefocaldata.R")
source("~/code/OrdRegMix/preparemodel.R")

std <- standardize(as.formula(paste0("SDB ~ ",leftside)),cov_dat,family=binomial)

Xf <- model.matrix(lm(formula=std$formula,data=std$data))[,-1]
Xr <- model.matrix( ~ 0 + FocalID,data=cov_dat) %>% as.matrix()
n <- used_obs[,length(Observation),by=.(FocalID,Year)]$V1
M <- length(n)

mvb <- stan_model("~/code/MultVarBinom/mvbinom.stan")
mvb_mm <- stan_model("~/code/MultVarBinom/mvbinom_mm.stan")
mvb_mm_f1 <- stan_model("~/code/MultVarBinom/mvbinom_mm_f1.stan")
# mvb_mm_f2 <- stan_model("~/code/MultVarBinom/mvbinom_mm_f2.stan")

standat <- list(Y=t(Y),N=nrow(Y),D=ncol(Y),M=length(n),n=n,D2=2^nrow(Y),X=t(Xf),P=ncol(Xf))

fit0 <- sampling(mvb,standat,
                 iter=534,warmup=200,thin=2,chains=3,cores=3)

standat <- c(standat,list(Z=t(Xr),R=nZ(Xr)))
# fit_mm_f2 <- sampling(mvb_mm_f2,standat,pars=c("f2_u_raw","f2_sigma_raw"),
#                       iter=700,warmup=200,thin=2,chains=4,include=F,cores=4)
fit_mm_f1 <- sampling(mvb_mm_f1,standat,pars=c("f1_u_raw"),
                      iter=534,warmup=200,thin=2,chains=3,include=F,cores=3)
fit_mm <- sampling(mvb_mm,standat,pars=c("f1_u_raw","f2_u_raw","f2_sigma_raw"),
                   iter=534,warmup=200,thin=2,chains=3,include=F,cores=3)
