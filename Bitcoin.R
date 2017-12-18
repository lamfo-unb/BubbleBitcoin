library(anytime)
library(tidyverse)
library(lubridate)
library(BayesianTools)

#https://www.coindesk.com/price/
btc<-read.csv("Data/BTC2016.csv")
btc$Date<-as.character(btc$Date)

#Log-likelihood
N <- nrow(btc)
Y <- log(btc$Close)
X <- seq(1,nrow(btc))

plot(X,Y,type="l")


loglik <- function(param) {
  ll <- 0
  A <- param[1]
  B <- param[2]
  beta<-param[3]
  tc <- param[4]
  C <- param[5]
  omega<- param[6]
  phi<- param[7]
  sigma<-param[8]
  for (i in 1:N) {
    mu<- A+(B*((tc-X[i])^beta))*(1+C*cos(omega*log(tc-X[i])+phi))
    ll <- ll + dnorm(Y[i], mu, sigma, log = TRUE)
  }
  return(ll)
}

maxTime<-max(X)
setup <- createBayesianSetup(loglik, prior = NULL,
                             lower = c(0, -100, 0.15, maxTime, -100, -4, 0, 0.01),
                             upper = c(100, -0.01, 0.51, 10000, 100, +8, 2*pi, +100))

out <- runMCMC(setup, sampler = "DE",
               settings = list(iterations = 100000,
                               burnin = 1500,
                               thin = 10))

gelmanDiagnostics(out)
summary(out)
plot(out)
correlationPlot(out)
marginalPlot(out)
MAP(out)
DIC(out)
WAIC(out) 
marginalLikelihood(out)


parms<-summary(out)
lambda <- 4737.526
kappa <- 36.639
beta <- 0.342
omega <- 5.581
phi<- 3.565
B0<- -27.852
B1<- -42.254
sigma<-53.494

ggplot(btc, aes(x = X, y = Y)) + geom_line()+
stat_function(fun=function(x) log(lambda)-(kappa/beta)*((B0*(lambda-x)^beta)+((B1*(lambda-x)^beta)*cos(omega*log(lambda-x)+phi))), geom="line") 




f <- ggplot(data.frame(x = X), aes(x))
f+stat_function(fun=function(x) log(lambda)-(kappa/beta)*((B0*(lambda-x)^beta)+((B1*(lambda-x)^beta)*cos(omega*log(lambda-x)+phi))), geom="line") 





density = function(param){
  lambda <- param[1]
  d1 <- dnorm(param[1], max(X),6, log =TRUE)
  
  
  
  kappa <- param[2]
  beta <- param[3]
  omega <- param[4]
  phi<- param[5]
  B0<- param[6]
  B1<- param[7]
  sigma<-param[8]
  
  d2 = dnorm(par[2], mean= 2, sd = 3, log =TRUE)
  return(d1 + d2)
}

# The sampling is optional but recommended because the MCMCs can generate automatic starting
# conditions if this is provided
sampler = function(n=1){
  d1 = runif(n, -2,6)
  d2 = rnorm(n, mean= 2, sd = 3)
  return(cbind(d1,d2))
}

prior <- createPrior(density = density, sampler = sampler, 
                     lower = c(-3,-3), upper = c(3,3), best = NULL)




setup <- createBayesianSetup(loglik, prior = prior,
                             lower = c(8000, -100, 0.15, 4, 0, -100, -100, 0.01),
                             upper = c(10000, +100, 0.5, 8, 2*pi, +100, +100, +100))

out <- runMCMC(setup, sampler = "DE",
               settings = list(iterations = 61500,
                               burnin = 1500,
                               thin = 5))

gelmanDiagnostics(out)
summary(out)
