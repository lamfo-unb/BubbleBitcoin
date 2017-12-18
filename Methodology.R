library(anytime)
library(tidyverse)
library(lubridate)
library(BayesianTools)
library(RcppDE)
library(foreach)
library(doParallel)

#Read the Data
coins.df<-readRDS("Data/Coins.RDS")

#Least Square
leastSquare<-function(parm)
{
  m<-parm[1]
  w<-parm[2]
  tc<-parm[3]
  phi<-parm[4]
  A<-parm[5]
  B<-parm[6]
  C<-parm[7]
  LS<-(Y0-(A+B*(tc-X0)^m*(1+C*cos(w*log(tc-X0)+phi))))^2
  #Constraints
  if(m>1 || m<0 || w<0||w>10|| phi<0||phi>2*pi||B>0 )
  {
    return(Inf)
  }
  else
  {
    return(sum(LS))  
  }
}

# Set seed
set.seed(5684)

#Coins
coin<-c("BTC","LTC","ETH","XRP","BCH","STR","XMR","XEM","ETC","DASH","STRAT","NXT","ZEC","EMC2","ARDR","LSK","DOGE","VTC","BTS","OMG","SC","DGB","FCT","MAID","FLDC","VRC","REP","ZRX","GNO","CVC","DCR","GNT","SYS","STEEM","GAME","LBC","NAV","BCN","AMP","CLAM","POT","STORJ","RIC","BURST","PPC","BLK","PASC","GAS","GRC","SBD","XCP","BTM","EXP","VIA","BTCD","RADS","FLO","NEOS","NMC","NXC","XBC","XVC","PINK","OMNI","BELA","BCY","HUC","XPM")
resultsBoot<-list()
for(i in 1:length(coin)){
  #Data
  Y <- log(as.numeric(unlist(coins.df[,coin[i]])))
  X <- seq(1,length(na.omit(Y)))
  
  #Maximum time
  minTime<-max(X)

  #Bootstrap
  ncl<-detectCores()
  cl <- makeCluster(ncl)
  registerDoParallel(cl)
  boot_b <- foreach(i=1:5000, .combine=rbind, .packages = c("RcppDE")) %dopar% {
    ind<- sample(length(X),length(X),replace=T)
    Y0 <- Y[ind]
    X0 <- X[ind]
    resultado2 <- DEoptim(leastSquare,
                          lower=c(0.33-0.18,6.36-1.56,minTime+2  ,0.01  ,0.01,-100,-0.999) ,
                          upper=c(0.33+0.18,6.36+1.56,10000,(2*pi-1e-10),100,-0.01, 0.999),control=list(NP=1000, itermax=400,trace=FALSE))
    parametros2<-c(resultado2$optim$bestmem)
    parametros2
  }
  
  #Stop clusters
  stopCluster(cl)

  #Save List
  resultsBoot[[i]]<-boot_b

}

save.image("Data/Bootstrap.RData")