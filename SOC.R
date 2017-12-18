library(anytime)
library(tidyverse)
library(lubridate)
library(BayesianTools)

#http://api.bitcoincharts.com/v1/csv/
btc<-read.csv("Data/btc24USD.csv")
colnames(btc)<-c("TIMESTAMP","PRICE","VOLUME")
btc$TIME<-anytime(btc$TIMESTAMP)

#Summarize by hour
btc2 <- btc %>% 
        mutate(hour = lubridate::floor_date(TIME , '5 hours')) %>% 
        group_by(hour) %>% 
        summarise(PRICE=mean(PRICE))
btc2<-na.omit(btc2)


#Log-likelihood
N <- nrow(btc2)
Y <- log(btc2$PRICE)
X <- as.numeric(difftime(btc2$hour,btc2$hour[1], units="hours"))

plot(X,Y,type="l")


loglik <- function(param) {
  ll <- 0
  lambda <- param[1]
  kappa <- param[2]
  beta <- param[3]
  omega <- param[4]
  phi<- param[5]
  B0<- param[6]
  B1<- param[7]
  sigma<-param[8]
  for (i in 1:N) {
    mu<- log(lambda)-(kappa/beta)*((B0*(lambda-X[i])^beta)+((B1*(lambda-X[i])^beta)*cos(omega*log(lambda-X[i])+phi)))
    ll <- ll + dnorm(Y[i],mu, sigma, log = TRUE)
  }
  return(ll)
}


setup <- createBayesianSetup(loglik, prior = NULL,
                             lower = c(8011, -100, 0.15, 4, 0, -100, -100, 0.01),
                             upper = c(10000, +100, 0.5, 8, 2*pi, +100, +100, +100))

out <- runMCMC(setup, sampler = "DE",
               settings = list(iterations = 100000,
                               burnin = 1500,
                               thin = 10))

gelmanDiagnostics(out)
summary(out)


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
