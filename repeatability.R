if (!exists("ffile")) ffile <- "mvbfits_rank.Rdat"
load(paste0("~/Dropbox/",ffile))
source('~/code/MultVarBinom/helperfuncs.R')

# corsd_mm <- calc_par_f(calc_cor_f,f1,f2,1) %>% 
#   apply(3,function(x) sd(x-rowMeans(x)))
# musd_mm <- calc_par_f(calc_marg_f,f1,f2,1) %>% 
#   apply(3,function(x) sd(x-rowMeans(x)))

pmm <- calc_p_mm(fit_mm)
rep_mm <- apply(pmm,3,calc_mi)

pf1 <- calc_p_f1(fit_mm_f1)
rep_f1 <- apply(pf1,3,calc_mi)

repdat <- data.table(rep=c(mean(rep_f1),mean(rep_mm)),
                     sd=c(sd(rep_f1),sd(rep_mm)),
                     model=c("Rates only","Rates & Cors"))
           # lb=c(quantile(rep_f1,0.025),quantile(rep_mm,0.025)),
           # ub=c(quantile(rep_f1,0.975),quantile(rep_mm,0.975)),

# ggplot(repdat,aes(y=Repeatability,x=Model,ymin=lb,ymax=ub)) + geom_pointrange()


# samps <- extract(fit_mm_f1,pars=c("f1_mu","f2"))
# p0 <- calc_p(t(samps$f1_mu),t(samps$f2)) %>% rowMeans()

# rr <- vector()
# for (i in 1:100) {
#   cat(i,"\r")
#   resamp <- sample.int(m,replace = T)
#   rmm <- apply(pmm[,resamp,],3,calc_repeat)
#   rf1 <- apply(pf1[,resamp,],3,calc_repeat)
#   rr[i] <- mean(rf1)/mean(rmm)
# }
