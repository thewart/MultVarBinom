simp <- function(d,musd,rhosd,n=1000,mu=0.0,rho=0.0) {
  require(parallel)
  require(mvtnorm)
  yset <- make_yset(d)
  ub <- lb <- yset
  ub[ub==-1] <- 0
  ub[ub==1] <- Inf
  lb[lb==-1] <- -Inf
  lb[lb==1] <- 0
  
  a <- rnorm(n*d,mu,musd) %>% matrix(nrow=d)
  b <- rnorm(n*d*(d-1)/2,rho,rhosd) %>% matrix(ncol=n)
  al <- plyr::alply(a,2)
  bl <- plyr::alply(b,2)
  
  return(mcmapply(simp_calcp,al,bl,MoreArgs = list(ub=ub,lb=lb), mc.cores = 6))
}

simp_calcp <- function(a,b,ub,lb) {
  d <- length(a)
  S <- diag(d)
  S[lower.tri(S)] <- S[upper.tri(S)] <- b
  S <- Matrix::nearPD(S)$mat %>% as.matrix()
  
  p <- vector(length=2^d)
  for (i in 1:length(p)) p[i] <- pmvnorm(lb[,i],ub[,i],a,sigma = S)
  return(p)
}

# for two dimensions -----
rhosd <- seq(0,0.4,0.05) %>% rep(3)
musd <- rep(c(0.025,0.125,0.26),each=length(rhosd)/3)
r2frac=vector()
corsd=vector()
mupsd=vector()
for (i in 1:length(musd)) {
  cat(i,"\r")
  p <- simp(2,musd[i],rhosd[i],mu=0,n=1000)
  p[p<0.001] <- 0.001
  p <- p/colSums(p)
  r2frac[i] <- (calc_repeat(p) - calc_repeat_marg(p))/calc_repeat(p)
  corsd[i] <- calc_cor(p) %>% sd()
  mupsd[i] <- calc_marg(p) %>% sd()
}
musdl <- rep(c(0.01,0.05,0.1),each=length(rhosd)/3)
cordat <- data.table(r2frac,corsd,musdl)
corplt <- ggplot(cordat,aes(y=r2frac,x=corsd,color=ordered(musdl))) + geom_line() + 
  coord_cartesian(ylim=c(0,1)) + xlab("Correlation SD") + ylab("% Repeatability") + 
  scale_color_discrete("Average SD")

#rho -> cor map for various d ----
d <- rep(2:6,each=5)
rhosd <- seq(0.01,0.21,0.05) %>% rep(5)
r2frac=vector()
corsd=vector()
mupsd=vector()
cormu <- vector()
for (i in 1:length(d)) {
  cat(i,"\r")
  p <- simp(d[i],0.125,rhosd[i],mu=0,n=1000)
  p[p<0.001] <- 0.001
  p <- p/colSums(p)
  r2frac[i] <- (calc_repeat(p) - calc_repeat_marg(p))/calc_repeat(p)
  corsd[i] <- calc_cor(p) %>% sd()
  mupsd[i] <- calc_marg(p) %>% sd()
  cormu[i] <- calc_cor(rowMeans(p) %>% matrix(nrow=nrow(p))) %>% mean()
}
pltdat <- data.table(rhosd,corsd,d)
ggplot(pltdat,aes(y=corsd,x=rhosd,color=ordered(d))) + geom_line() + geom_point()
pltlm <- pltdat[,.(lm(corsd ~ rhosd)$coef[1],lm(corsd ~ rhosd)$coef[2]),by=d]
rhosd <- sapply(c(0.025,0.05,0.1),function(x) pltlm[,(x-V1)/V2]) %>% as.vector()

# across d for realsies ----
#rhosd <- c(c(0.02,0.024,0.025) #get rhosd from above
d <- rep(2:6,3)
r2frac <- vector()
corsd <- vector()
mupsd <- vector()
cormu <- vector()
for (i in 1:length(d)) {
  cat(i,"\r")
  p <- simp(d[i],0.125,rhosd[i],mu=0)
  p[p<0.001] <- 0.001
  p <- p/colSums(p)
  cormu[i] <- calc_cor(rowMeans(p)%>%matrix(nrow=nrow(p))) %>% mean()
  r2frac[i] <- (calc_repeat(p) - calc_repeat_marg(p))/calc_repeat(p)
  corsd[i] <- calc_cor(p) %>% sd()
  mupsd[i] <- calc_marg(p) %>% sd()
}
corsdl <- rep(c(0.01,0.05,0.1),each=length(rhosd)/3)
ddat <- data.table(r2frac=c(0,0,0,r2frac),d=c(1,1,1,d),corsdl=c(0.01,0.05,0.1,corsdl))
#color <- ggplot_build(corplt)$data[[1]]$colour %>% unique
dimplt <- ggplot(ddat,aes(y=r2frac,x=d,color=ordered(corsdl))) + geom_line() + geom_point() +
  scale_x_continuous("Number of behaviors",breaks=1:6,labels=as.character(1:6),minor_breaks = NULL) +
  ylab("% Repeatability") + coord_cartesian(ylim=c(0,1)) + scale_color_manual("Correlation SD",values=c("#E69F00", "#56B4E9", "#009E73"))


