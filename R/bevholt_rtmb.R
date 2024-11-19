library(RTMB)

source("../data/recdata.R")

objfun <- function(par, data)
{
  Rmax <- exp(par[["logRmax"]])
  S50 <- exp(par[["logS50"]])
  sigma <- exp(par[["logSigma"]])

  R <- data$R
  S <- data$S

  Rhat <- Rmax * S / (S + S50)
  nll <- -sum(dnorm(log(R), log(Rhat), sigma, TRUE))
  ADREPORT(log(Rhat))
  nll
}

bevholt <- function(f, data) function(par) f(par, data)

par <- list(logRmax=0, logS50=0, logSigma=0)
model <- MakeADFun(bevholt(objfun, recdata), par, silent=TRUE)
fit <- nlminb(model$par, model$fn, model$gr)
rep <- summary(sdreport(model))

hessian <- model$he()
SE <- sqrt(diag(solve(hessian)))
cor <- round(solve(hessian) / (SE %o% SE), 3)

rep
cor
fit
