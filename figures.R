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
mu_std <- apply(yhat,c(1,3),sd)
ycor <- calc_par(calc_cor,f1,f2)
cor_std <- apply(ycor,c(1,3),sd)

mudat_mm <- make_pltdat(mu_std,label = "mm")
cordat_mm <- make_pltdat(cor_std,label = "mm")



samps <- extract(fit_mm_f1,pars=c("f1_mu","f1_u","f2"))
f1 <- array(dim=dim(samps$f1_u))
for (i in 1:n) f1[i,,] <- samps$f1_mu[i,] + samps$f1_u[i,,]
f2 <- array(samps$f2,dim=c(dim(samps$f2),m))

yhat <- calc_par(calc_marg,f1,f2,cores=4)
mu_std <- apply(yhat,c(1,3),sd)
ycor <- calc_par(calc_cor,f1,f2)
cor_std <- apply(ycor,c(1,3),sd)

mudat_f1 <- make_pltdat(mu_std,label="f1")
cordat_f1 <- make_pltdat(cor_std,label = "f1")




samps <- extract(fit_mm_f2,pars=c("f1_mu","f2_mu","f2_u"))
f2 <- array(dim=dim(samps$f2_u))
for (i in 1:n) f2[i,,] <- samps$f2_mu[i,] + samps$f2_u[i,,]
f1 <- array(samps$f1_mu,dim=c(dim(samps$f1_mu),m))
yhat <- calc_par(calc_marg,f1,f2,cores=2)
mu_std <- apply(yhat,c(1,3),sd)
ycor <- calc_par(calc_cor,f1,f2)
cor_std <- apply(ycor,c(1,3),sd)

mudat_f2 <- make_pltdat(mu_std,label="f2")
cordat_f2 <- make_pltdat(cor_std,label = "f2")

bigdat <- rbindlist(list(mudat_mm,mudat_f1,mudat_f2))
bigdat[,behav:=behaviors]

mydodge <- position_dodge(width=0.5)
ggplot(bigdat,aes(y=std,x=behav,color=`label`),position) + geom_point(size=2,position=mydodge) + 
  # geom_errorbar(aes(ymin=lbi,ymax=ubi),width=0,size=1,position=mydodge) +
  geom_errorbar(aes(ymin=lbo,ymax=ubo),width=0,position=mydodge) + coord_flip()
