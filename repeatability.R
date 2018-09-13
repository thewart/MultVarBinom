X <- model.matrix( ~ 0 + FocalID,data=used_obs) %>% as.matrix() %>% t()
lp0 <- extract(fit0,"lp")$lp %>% colMeans()
lp_mm <- extract(fit_mm,"lp")$lp %>% colMeans()
lp_f1 <- extract(fit_mm_f1,"lp")$lp %>% colMeans()

Rmm <- -(X %*% (lp_mm-lp0)) / (X %*% lp0)
Rf1 <- -(X %*% (lp_f1-lp0)) / (X %*% lp0)
Rdiff <- (Rmm-Rf1)/Rmm
Rdiff[Rdiff<0] <- 0 #hack for correcting monte-carlo error
