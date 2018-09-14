X <- model.matrix( ~ 0 + FocalID,data=used_obs) %>% as.matrix() %>% t()
lp0 <- extract(fit0,"lp")$lp %>% colMeans()
lp_mm <- extract(fit_mm,"lp")$lp %>% colMeans()
lp_f1 <- extract(fit_mm_f1,"lp")$lp %>% colMeans()
lp_f2 <- extract(fit_mm_f2,"lp")$lp %>% colMeans()


Rmm <- (X %*% (lp_mm-lp0)) / rowSums(X)
Rf1 <- (X %*% (lp_f1-lp0)) / rowSums(X)
R0 <- (X %*% -lp0) / rowSums(X)
Rdiff <- (X %*% (lp_mm-lp_f1)) / rowSums(X)
RDR <- Rdiff/Rmm 
RDR[RDR<0] <- 0 #hack for correcting monte-carlo error
