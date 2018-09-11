samps <- extract(fit_mm,pars=c("f1","f2"))
n <- dim(samps$f1)[1]
d <- dim(samps$f1)[2]
m <- dim(samps$f1)[3]
yhat <- calc_par(calc_marg,samps$f1,samps$f2,cores = 2)
mu_std <- apply(yhat,c(1,3),sd)
ycor <- calc_par(calc_cor,samps$f1,samps$f2)
cor_std <- apply(ycor,c(1,3),sd)

mudat_mm <- make_pltdat(mu_std,label = "mm")
cordat_mm <- make_pltdat(cor_std,label = "mm")

samps <- extract(fit_mm_f1,pars=c("f1","f2"))
samps$f2 <- array(samps$f2,dim=c(dim(samps$f2),m))
yhat <- calc_par(calc_marg,samps$f1,samps$f2,cores=2)
mu_std <- apply(yhat,c(1,3),sd)
ycor <- calc_par(calc_cor,samps$f1,samps$f2)
cor_std <- apply(ycor,c(1,3),sd)

mudat_f1 <- make_pltdat(mu_std,label="f1")
cordat_f1 <- make_pltdat(cor_std,label = "f1")

samps <- extract(fit_mm_f2,pars=c("f1","f2"))
yhat <- calc_par(calc_marg,samps$f1,samps$f2,cores=2)
mu_std <- apply(yhat,c(1,3),sd)
ycor <- calc_par(calc_cor,samps$f1,samps$f2)
cor_std <- apply(ycor,c(1,3),sd)

mudat_f2 <- make_pltdat(mu_std,label="f2")
cordat_f2 <- make_pltdat(cor_std,label = "f2")

bigdat <- rbindlist(list(cordat_mm,cordat_f1,cordat_f2))
bigdat[,behav:=behaviors]

mydodge <- position_dodge(width=0.5)
ggplot(bigdat,aes(y=std,x=rep(factor(1:45),3),color=`label`),position) + geom_point(size=2,position=mydodge) + 
  # geom_errorbar(aes(ymin=lbi,ymax=ubi),width=0,size=1,position=mydodge) +
  geom_errorbar(aes(ymin=lbo,ymax=ubo),width=0,position=mydodge) + coord_flip()
