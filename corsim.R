d <- 5
n <- 100
mu <- 0.0
musd <- 0.1
rho <- 0.0
rhosd <- 0.1

yset <- make_yset(d)
ub <- lb <- yset
ub[ub==-1] <- 0
ub[ub==1] <- Inf
lb[lb==-1] <- -Inf
lb[lb==1] <- 0

ubl <- plyr::alply(ub,2)
lbl <- plyr::alply(lb,2)

a <- rnorm(n*d,mu,musd) %>% matrix(nrow=d)
b <- rnorm(n*d*(d-1)/2,rho,rhosd) %>% matrix(ncol=n)

al <- plyr::alply(a,2)
bl <- plyr::alply(b,2)

# p <- matrix(nrow=2^d,ncol=n)
# for (i in 1:n) {
#   cat(i,"\r")
#   S <- diag(d)
#   S[lower.tri(S)] <- S[upper.tri(S)] <- b[,i]
#   S <- Matrix::nearPD(S)$mat %>% as.matrix()
#   
#    for (j in 1:nrow(p)) p[j,i] <- pmvnorm(lb[,j],ub[,j],a[,i],sigma = S)
#   #p[,i] <- mcmapply(pmvnorm,lbl,ubl,MoreArgs = list(mean=a[,i],sigma=S),mc.cores=4)
# }

p <- 

simp <- function(a,b,ub,lb) {
  d <- length(a)
  S <- diag(d)
  S[lower.tri(S)] <- S[upper.tri(S)] <- b
  S <- Matrix::nearPD(S)$mat %>% as.matrix()
  
  p <- vector(length=2^d)
  for (i in 1:length(p)) p[i] <- pmvnorm(lb[,i],ub[,i],a,sigma = S)
  return(p)
}
