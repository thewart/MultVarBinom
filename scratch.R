if (!exists("ffile")) ffile <- "mvbfits_rank.Rdat"
load(paste0("~/Dropbox/",ffile))
source("~/code/MultVarBinom/helperfuncs.R")

behnames <- c("SDB","Feed","Travel","PassiveContact","Agg_give","Avoid","Submit",
  "Groom_high","Groom_low","Approach_high","Approach_low")
pmm <- calc_p_mm(fit_mm)
mu <- calc_par(calc_marg,pmm)
cors <- calc_par(calc_cor,pmm)

foo <- apply(cors,c(1,3),mean)
pltdt <- mat2gg(foo)
pltdt$behav <- interact_string(behnames)
ggplot(pltdt,aes(y=reorder(behav,mu),x=mu,xmin=lb,xmax=ub)) + 
  geom_point() + geom_errorbarh(height=0) + ylab(NULL) + xlab("Correlation") + theme_light()

juh <- diag(length(behnames))*0.5
juh[lower.tri(juh)] <- pltdt$mu
juh <- juh + t(juh)
rownames(juh) <- colnames(juh) <- behnames
mat2gg <- function(x,p=c(0.025,0.975)) {
  data.table(mu=rowMeans(x),
             sd=apply(x,1,sd),
             lb=apply(x,1,quantile,probs=p[1]),
             ub=apply(x,1,quantile,probs=p[2]))
}

kldist <- function(p,vec=F) {
  n <- ncol(p)
  kldmat <- array(dim = c(n,n))
  for (i in 1:n) {
    for (j in i:n) {
      ld <- log(p[,i])-log(p[,j])
      kldmat[i,j] <- kldmat[j,i] <- sum(p[,i]*ld - p[,j]*ld)
    }
  }
  if (vec) kldmat <- kldmat[lower.tri(kldmat)]
  return(kldmat)
}

foo <- apply(pmm,c(3),function(x) mean(kldist(x,T)))
juh <- apply(pmm,c(3),function(x) mean(kldist(match_marg(x),T)))