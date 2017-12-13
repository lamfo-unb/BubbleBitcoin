library(anytime)
library(tidyverse)
library(lubridate)
library(BayesianTools)
library(DEoptim)
library(foreach)
library(doParallel)

#https://www.coindesk.com/price/
btc<-read.csv("Data/BTC2016.csv")
btc$Date<-as.character(btc$Date)
btc<- na.omit(btc)

#Log-likelihood
N <- nrow(btc)
Y <- log(btc$Close)
X <- seq(1,nrow(btc))

plot(X,Y,type="l")

leastSquare<-function(parm)
{
  m<-parm[1]
  w<-parm[2]
  tc<-parm[3]
  phi<-parm[4]
  A<-parm[5]
  B<-parm[6]
  C<-parm[7]
  LS<-(Y-(A+B*(tc-X)^m*(1+C*cos(w*log(tc-X)+phi))))^2
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

set.seed(1234)
minTime<-max(X)
resultado2 <- DEoptim(leastSquare,
                      lower=c(0.33-0.18,6.36-1.56,minTime+2  ,0.01  ,0.01,-100,-0.999) ,
                      upper=c(0.33+0.18,6.36+1.56,10000,(2*pi-1e-10),100,-0.01, 0.999),control=list(NP=1000, itermax=400,trace=FALSE))
parametros2<-c(resultado2$optim$bestmem)


#Bootstrap
btc$X<-seq(1,nrow(btc))
ncl<-detectCores()
cl <- makeCluster(ncl)
registerDoParallel(cl)
boot_b <- foreach(i=1:10, .combine=rbind, .packages = c("DEoptim")) %dopar% {
  bdados<- btc[sample(nrow(btc),nrow(btc),replace=T),]
  N <- nrow(bdados)
  Y <- log(bdados$Close)
  X <- bdados$X
  resultado2 <- DEoptim(leastSquare,
                        lower=c(0.33-0.18,6.36-1.56,minTime+2  ,0.01  ,0.01,-100,-0.999) ,
                        upper=c(0.33+0.18,6.36+1.56,10000,(2*pi-1e-10),100,-0.01, 0.999),control=list(NP=1000, itermax=400,trace=FALSE))
  parametros2<-c(resultado2$optim$bestmem)
  parametros2
}

#Stop clusters
stopCluster(cl)


#Bootstrap estimates
parms.means<-colMeans(boot_b)
parms.sd <-apply(boot_b,1,sd)

m<-parms.means[1]
w<-parms.means[2]
tc<-parms.means[3]
phi<-parms.means[4]
A<-parms.means[5]
B<-parms.means[6]
C<-parms.means[7]

time<-seq(ymd('2016-01-01'), length.out = 2000 , by = '1 day')
data$Date<-time

data<-data.frame(X=c(X,rep(NA,2000-length(X))),Y=c(Y,rep(NA,2000-length(Y))),Date=time)

ggplot(data, aes(x = X, y = Y)) +  geom_line() +
  stat_function(fun=function(x) A+B*(tc-x)^m*(1+C*cos(w*log(tc-x)+phi)), geom="line", color="red",  xlim = c(0,2000))+
  geom_vline(xintercept=1750)+
  labs(x="Time",y="Price") +
  annotate("text", x = 1750, y=13, label = "15 de Outubro de 2020")

ggplot(na.omit(data), aes(x = Date, y = Y)) +  geom_line() +labs(x="Time",y="Price") 


