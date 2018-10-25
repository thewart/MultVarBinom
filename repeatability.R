#load("~/Dropbox/mvbfits.Rdat")
library(parallel)
samps <- extract(fit_mm,pars=c("f1_mu","f1_u","f2_mu","f2_u"))
n <- dim(samps$f1_u)[1]
d <- dim(samps$f1_u)[2]
m <- dim(samps$f1_u)[3]
f1 <- array(dim=dim(samps$f1_u))
f2 <- array(dim=dim(samps$f2_u))
for (i in 1:n) {
  f1[i,,] <- samps$f1_mu[i,] + samps$f1_u[i,,]
  f2[i,,] <- samps$f2_mu[i,] + samps$f2_u[i,,]
}

corsd_mm <- calc_par_f(calc_cor_f,f1,f2,6) %>% 
  apply(3,function(x) sd(x-rowMeans(x)))
musd_mm <- calc_par_f(calc_marg_f,f1,f2,6) %>% 
  apply(3,function(x) sd(x-rowMeans(x)))

pmm <- calc_par_f(calc_p,f1,f2,6)
rep_mm <- apply(pmm,3,calc_repeat)

samps <- extract(fit_mm_f1,pars=c("f1_mu","f1_u","f2"))
f1 <- array(dim=dim(samps$f1_u))
for (i in 1:n) f1[i,,] <- samps$f1_mu[i,] + samps$f1_u[i,,]
f2 <- array(samps$f2,dim=c(dim(samps$f2),m))

pf1 <- calc_par_f(calc_p,f1,f2,6)
rep_f1 <- apply(pf1,3,calc_repeat)

repdat <- data.table(Repeatability=c(mean(rep_f1),mean(rep_mm)),
           lb=c(quantile(rep_f1,0.025),quantile(rep_mm,0.025)),
           ub=c(quantile(rep_f1,0.975),quantile(rep_mm,0.975)),
           Model=c("Rates","Rates+Interactions"))
ggplot(repdat,aes(y=Repeatability,x=Model,ymin=lb,ymax=ub)) + geom_pointrange() + 
  ylim(c(0,0.0325))

# samps <- extract(fit_mm_f1,pars=c("f1_mu","f2"))
# p0 <- calc_p(t(samps$f1_mu),t(samps$f2)) %>% rowMeans()

rr <- vector()
for (i in 1:100) {
  cat(i,"\r")
  resamp <- sample.int(94,replace = T)
  rmm <- apply(pmm[,resamp,],3,calc_repeat)
  rf1 <- apply(pf1[,resamp,],3,calc_repeat)
  rr[i] <- mean(rf1)/mean(rmm)
}
