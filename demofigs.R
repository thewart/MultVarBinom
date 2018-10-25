library(RColorBrewer)
library(ggExtra)
library(patchwork)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0, 1))
scnull <- scale_fill_gradientn(colours = myPalette(100),limits=c(0,1),guide=F)
nomarg <- theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
make_marg_plt <- function(groom,bite,scale=F) {
  gd <- data.table(Groom=ordered(c(0,1)),Bite=ordered(c(0,0)),Rate=c(1-groom,groom))
  bd <- data.table(Groom=ordered(c(0,0)),Bite=ordered(c(0,1)),Rate=c(1-bite,bite))
  
  gplt <- ggplot(gd,aes(y=Groom,x=Bite,fill=Rate)) + geom_tile() + 
    theme_minimal() + coord_equal() + scnull + nomarg + 
    scale_x_discrete(NULL,labels=NULL) + removeGrid() 
  bplt <- ggplot(bd,aes(y=Bite,x=Groom,fill=Rate)) + geom_tile() +
    theme_minimal() + coord_equal() + nomarg +
    scale_x_discrete(NULL,labels=NULL) + removeGrid()
  if (scale) { bplt <- bplt + sc
  } else { bplt <- bplt + scnull}
  
  return(list(gplt,bplt))
}

make_joint_plt <- function(groom,bite,cor,scale=F) {
  #check bounds
  lb <- max( -sqrt( groom*bite/((1-groom)*(1-bite)) ),
             -sqrt( (1-groom)*(1-bite)/(groom*bite) ))
  ub <- min( sqrt( groom*(1-bite)/(bite*(1-groom)) ), 
             sqrt( bite*(1-groom)/(groom*(1-bite)) ))
  if (cor<lb) cor <- lb
  if (cor>ub) cor <- ub
  
  p11 <- cor*sqrt(groom*(1-groom)*(1-bite)*bite) + groom*bite
  p10 <- groom - p11
  p01 <- bite - p11
  p00 <- 1 - p11 - p10 - p01
  
  dat <- data.table(Groom=ordered(c(0,1,0,1)),
                    Bite=ordered(c(0,0,1,1)),
                    Rate=c(p00,p10,p01,p11))
  
  jplt <- ggplot(dat,aes(y=Groom,x=Bite,fill=Rate)) + geom_tile() + 
    theme_minimal() + coord_equal() + removeGrid() + nomarg
  
  if (scale) { jplt <- jplt + sc
  } else { jplt <- jplt + scnull}
  
  return(jplt)
}

pop <- make_marg_plt(0.3,0.2,T)
plot_grid(pop[[1]],pop[[2]],rel_widths = c(1,1.5))

ind1 <- make_marg_plt(0.32,0.18)
ind2 <- make_marg_plt(0.1,0.1)
ind3 <- make_marg_plt(0.5,0.3)
ind4 <- make_marg_plt(0.28,0.22)
plot_grid(ind1[[1]],ind1[[2]],NULL,ind2[[1]],ind2[[2]],
          ind3[[1]],ind3[[2]],NULL,ind4[[1]],ind4[[2]],
          nrow=2,rel_widths = c(1,1,0.5,1,1),
          labels = c("A","","","B","","C","","","D"))


inda <- pop <- make_marg_plt(0.3,0.2,T)
ind1 <- make_joint_plt(0.3,0.2,0.0)
ind2 <- make_joint_plt(0.3,0.2,0.7)
ind3 <- make_joint_plt(0.3,0.2,-0.3)
plot_grid(ind1,ind2,ind3,nrow=1,labels=c("A1","A2","A3"))
