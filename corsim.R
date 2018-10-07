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

# for two dimensions
rhosd <- c(seq(0,0.75,0.05),seq(0.9,1.2,0.15)) %>% rep(3)
musd <- rep(c(0.175,0.35,0.75),each=length(rhosd)/3)
r2frac=vector()
corsd=vector()
mupsd=vector()
for (i in 1:length(musd)) {
  cat(i,"\r")
  p <- simp(2,musd[i],rhosd[i],mu=-0.84,n=1000)
  p[p<0.001] <- 0.001
  p <- p/colSums(p)
  r2frac[i] <- (calc_repeat(p) - calc_repeat_marg(p))/calc_repeat(p)
  corsd[i] <- calc_cor(p) %>% sd()
  mupsd[i] <- calc_marg(p) %>% sd()
}
musdl <- rep(c(0.05,0.1,0.2),each=length(rhosd)/3)
cordat <- data.table(r2frac,corsd,musdl)
corplt <- ggplot(cordat[corsd<=0.4],aes(y=r2frac,x=corsd,color=ordered(musdl))) + geom_line() + 
  coord_cartesian(ylim=c(0,1)) + xlab("Correlation SD") + ylab("% Repeatability") + 
  scale_color_discrete("Average SD")

#identify 0.1 point for various d ----
# d <- rep(2:6,each=6)
# rhosd <- seq(0.2,0.3,0.02) %>% rep(5)
# r2frac=vector()
# corsd=vector()
# mupsd=vector()
# for (i in 1:length(d)) {
#   cat(i,"\r")
#   p <- simp(d[i],0.35,rhosd[i],mu=-0.84,n=1000)
#   p[p<0.001] <- 0.001
#   p <- p/colSums(p)
#   r2frac[i] <- (calc_repeat(p) - calc_repeat_marg(p))/calc_repeat(p)
#   corsd[i] <- calc_cor(p) %>% sd()
#   mupsd[i] <- calc_marg(p) %>% sd()
# }
# pltdat <- data.table(rhosd,corsd,d)
# ggplot(pltdat,aes(y=corsd,x=rhosd,color=ordered(d))) + geom_line() + geom_point()


# across d for realsies ----
d <- 2:6
rhosd <- c(0.2,0.21,0.225,0.255,0.275)
r2frac <- vector()
corsd <- vector()
mupsd <- vector()
cormu <- vector()
for (i in 1:length(d)) {
  cat(i,"\r")
  p <- simp(d[i],0.35,rhosd[i],mu=-0.84)
  p[p<0.001] <- 0.001
  p <- p/colSums(p)
  cormu[i] <- calc_cor(rowMeans(p)%>%matrix(nrow=nrow(p))) %>% mean()
  r2frac[i] <- (calc_repeat(p) - calc_repeat_marg(p))/calc_repeat(p)
  corsd[i] <- calc_cor(p) %>% sd()
  mupsd[i] <- calc_marg(p) %>% sd()
}
ddat <- data.table(r2frac=c(0,r2frac),d=c(1,d))
color <- ggplot_build(corplt)$data[[1]]$colour %>% unique
dimplt <- ggplot(ddat,aes(y=r2frac,x=d)) + geom_line(color=color[2]) + coord_cartesian(ylim=c(0,1)) +
  scale_x_continuous("Number of behaviors",breaks=1:6,labels=as.character(1:6),minor_breaks = NULL) + ylab("% Repeatability")


