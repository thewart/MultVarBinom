if (!exists("ffile")) ffile <- "mvbfits_rank11.Rdat"
load(paste0("~/Dropbox/",ffile))

#####  waic
library(loo)
waic_mm <- extract_log_lik(fit_mm,"lp") %>% waic()
waic_f1 <- extract_log_lik(fit_mm_f1,"lp") %>% waic()
waic0 <- extract_log_lik(fit0,"lp") %>% waic()

waicdat <- rbind(c(0,0),compare(waic0,waic_f1),compare(waic0,waic_mm))%>% as.data.table()
waicdat$Model <- c("Nothing","Rates only","Rates & Cors")

# ggplot(waicdat[Model!="None"],aes(y=elpd_diff,x=Model)) + geom_line() + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se))
# ggplot(waicdat,aes(y=elpd_diff,x=Model)) + geom_line() + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se))

##### mutual information
source('~/code/MultVarBinom/helperfuncs.R')

pmm <- calc_p_mm(fit_mm)
rep_mm <- apply(pmm,3,calc_mi,normalize=F)
rep_mm_marg <- apply(pmm,3,calc_mi_marg,normalize=F)
pf1 <- calc_p_f1(fit_mm_f1)
rep_f1 <- apply(pf1,3,calc_mi,normalize=F)
rep_f1_marg <- apply(pf1,3,calc_mi_marg,normalize=F)

midiff_qeb <- (rep_mm-rep_mm_marg)/rep_mm_marg

