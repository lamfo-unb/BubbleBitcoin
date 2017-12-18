library(BayesianTools)

set.seed(123)

data <- read.csv("https://raw.githubusercontent.com/ito4303/naro_toukei/master/example4.csv")

N <- nrow(data)
Y <- data$Num
X <- data$Light

loglik <- function(param) {
  # param[1]: inclusion probability
  # param[2]: intercept of Poisson regression
  # param[3]: slope of Poisson regression
  ll <- 0
  lambda <- exp(param[2] + param[3] * X)
  for (i in 1:N) {
    if (Y[i] > 0) {
      ll <- ll +
        dbinom(1, 1, param[1], log = TRUE) +
        dpois(Y[i], lambda[i], log = TRUE)
    } else {
      p1 <- dbinom(0, 1, param[1])
      p2 <- dbinom(1, 1, param[1]) * dpois(0, lambda[i])
      ll <- ll + log(p1 + p2)
    }
  }
  return(ll)
}

setup <- createBayesianSetup(loglik, prior = NULL,
                             lower = c(0, -100, -100),
                             upper = c(1, 100, 100))
out <- runMCMC(setup, sampler = "DEzs",
               settings = list(iterations = 61500,
                               burnin = 1500,
                               thin = 5))


gelmanDiagnostics(out)
summary(out)

library(BayesianTools)
settings <- list(iterations = iter,  message = FALSE)
out <- runMCMC(bayesianSetup = bayesianSetup, sampler = "DE", settings = settings)
plot(out) 