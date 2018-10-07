#load("~/Dropbox/mvbfits.Rdat")
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
# pmm <- calc_par_f(calc_p,f1,f2,6) %>% apply(c(1,2),mean)
# cor_mm <- calc_cor(pmm)
# corsd_mm <- (cor_mm - rowMeans(cor_mm) ) %>% sd
# mu_mm <- calc_marg(pmm)
# musd_mm <- (mu_mm - rowMeans(mu_mm) ) %>% sd

corsd_mm <- calc_par_f(calc_cor_f,f1,f2,6) %>% 
  apply(3,function(x) sd(x-rowMeans(x)))
musd_mm <- calc_par_f(calc_marg_f,f1,f2,6) %>% 
  apply(3,function(x) sd(x-rowMeans(x)))

rep_mm <- calc_par_f(calc_p,f1,f2,6)


samps <- extract(fit_mm_f1,pars=c("f1_mu","f1_u","f2"))
f1 <- array(dim=dim(samps$f1_u))
for (i in 1:n) f1[i,,] <- samps$f1_mu[i,] + samps$f1_u[i,,]
f2 <- array(samps$f2,dim=c(dim(samps$f2),m))
pf1 <- calc_par(calc_p,f1,f2,6) %>% apply(c(1,2),mean)

pf1 <- calc_par(calc_p,f1,f2,6)
foo <- apply(pf1,3,calc_repeat)


samps <- extract(fit_mm_f1,pars=c("f1_mu","f2"))
p0 <- calc_p(t(samps$f1_mu),t(samps$f2)) %>% rowMeans()


rr <- vector()
for (i in 1:1000) {
  cat(i,"\r")
  resamp <- sample.int(94,replace = T)
  p0 <- rowMeans(pmm[,resamp])/2 + rowMeans(pf1[,resamp])/2
  rr[i] <- (calc_repeat(pmm[,resamp],p0)-calc_repeat(pf1[,resamp],p0))/calc_repeat(pmm[,resamp],p0)
}

pmm <- calc_par(calc_p,f1,f2,6)
juh <- apply(pmm,1,calc_repeat)
