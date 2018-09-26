crossent <- function(q,p,mu,yset) {
  return( -sum(p*log(q)) )
}
crossent_grad <- function(q,p,mu,yset) {
  return( -p/q )
}

crossent_const <- function(q,p,mu,yset) {
  return( as.vector(yset %*% q) - mu )
}

crossent_const_grad <- function(q,p,mu,yset) {
  return(yset)
}

mu <- c(0.3,0.3)
p <- calcp(0.2,0.4,0.3)
q <- calcp(0.3,0.3,0.0)

X <- log2(length(p)) %>% make_yset
X[X==-1] <- 0
X <- rbind(1,X)

nloptr(x0=q,
       eval_f = crossent,
       eval_grad_f = crossent_grad,
       lb = rep(0,length(p)),
       eval_g_eq = crossent_const,
       eval_jac_g_eq = crossent_const_grad,
       opts = list(algorithm="NLOPT_LD_SLSQP"),
       p = p,
       mu = c(1,mu),
       yset = X)

local_opts <- list( "algorithm" = "NLOPT_LD_LBFGS",
                    "xtol_rel" = 1.0e-4 )

opts <- list( "algorithm" = "NLOPT_LD_AUGLAG",
              "xtol_rel" = 1.0e-4,
              "maxeval" = 1000,
              "local_opts" = local_opts )
