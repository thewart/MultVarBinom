fit <- fit_mm
samps <- extract(fit,pars=c("f1","f2"))
n <- dim(samps$f1)[1]
d <- dim(samps$f1)[2]
m <- dim(samps$f1)[3]
yset <- make_yset(dim(samps$f1),[2])
xset <- make_xset(yset)

mu_std <- matrix(nrow=d,ncol=n)
cor_std <- matrix(nrow=d*(d-1)/2,ncol=n)
for (i in 1:n) {
  cat(i,"\r")
  p <- calc_p(samps$f1[i,,],samps$f2[i,,],yset,xset)
  yhat <- calc_marg(samps$f1[i,,],samps$f2[i,,],yset,p)
  ycor <- calc_cor(samps$f1[i,,],samps$f2[i,,],yset,p)
  mu_std[,i] <- apply(yhat,1,sd)
  cor_std[,i] <- apply(ycor,1,sd)
}

mudat <- data.table(std=apply(mu_std,1,median),
           ubi=apply(mu_std,1,quantile,0.1),ubo=apply(mu_std,1,quantile,0.025),
           lbi=apply(mu_std,1,quantile,0.9),lbo=apply(mu_std,1,quantile,0.975))
ggplot(mudat,aes(x=std,y=behaviors)) + geom_point(size=2) + 
  geom_errorbarh(aes(xmin=lbi,xmax=ubi),height=0,size=1.5) + geom_errorbarh(aes(xmin=lbo,xmax=ubo),height=0)

cordat <- data.table(std=apply(cor_std,1,median),
                    ubi=apply(cor_std,1,quantile,0.1),ubo=apply(cor_std,1,quantile,0.025),
                    lbi=apply(cor_std,1,quantile,0.9),lbo=apply(cor_std,1,quantile,0.975))
ggplot(cordat,aes(x=std,y=factor(1:45))) + geom_point(size=2) + 
  geom_errorbarh(aes(xmin=lbi,xmax=ubi),height=0,size=1.5) + geom_errorbarh(aes(xmin=lbo,xmax=ubo),height=0)
