#load("~/Dropbox/mvbfits.Rdat")
library(loo)
lp_mm <- extract_log_lik(fit_mm,"lp")
lp_f1 <- extract_log_lik(fit_mm_f1,"lp")
lp0 <- extract_log_lik(fit0,"lp")

waic_mm <- waic(lp_mm)
waic_f1 <- waic(lp_f1)
waic0 <- waic(lp0)

waicdat <- rbind(c(0,0),compare(waic0,waic_f1),compare(waic0,waic_mm))%>% as.data.table()
waicdat$Model <- c("None","Rates","Rates+Interactions")

ggplot(waicdat[Model!="None"],aes(y=elpd_diff,x=Model)) + geom_line() + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se))
ggplot(waicdat,aes(y=elpd_diff,x=Model)) + geom_line() + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se))
