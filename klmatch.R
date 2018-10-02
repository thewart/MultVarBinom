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

yhat <- calc_par(calc_marg,f1,f2,cores = 4)
phat <- calc_par(calc_p,f1,f2,cores = 4)

pave <- apply(phat,c(1,3),mean)

p <- rowMeans(phat[,,1])
i <- 10
q <- phat[,i,1]
mu <- yhat[,i,1]
X <- log2(length(p)) %>% make_yset
X[X==-1] <- 0
X <- rbind(1,X)

nloptr(x0=q,
       eval_f = crossent,
       eval_grad_f = crossent_grad,
       lb = rep(0,length(p)),
       eval_g_eq = crossent_const,
       eval_jac_g_eq = crossent_const_grad,
       opts = list(algorithm="NLOPT_LD_SLSQP",maxeval=1e6),
       p = p,
       mu = c(1,mu),
       yset = X)
