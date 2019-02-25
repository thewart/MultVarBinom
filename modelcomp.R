if (!exists("ffile")) ffile <- "mvbfits_rank.Rdat"
load(paste0("~/Dropbox/",ffile))
library(loo)
waic_mm <- extract_log_lik(fit_mm,"lp") %>% waic()
waic_f1 <- extract_log_lik(fit_mm_f1,"lp") %>% waic()
waic0 <- extract_log_lik(fit0,"lp") %>% waic()

waicdat <- rbind(c(0,0),compare(waic0,waic_f1),compare(waic0,waic_mm))%>% as.data.table()
waicdat$model <- c("Nothing","Rates only","Rates & Cors")

# ggplot(waicdat[Model!="None"],aes(y=elpd_diff,x=Model)) + geom_line() + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se))
# ggplot(waicdat,aes(y=elpd_diff,x=Model)) + geom_line() + geom_pointrange(aes(ymin=elpd_diff-se,ymax=elpd_diff+se))
