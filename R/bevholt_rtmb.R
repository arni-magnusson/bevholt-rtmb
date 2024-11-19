library(RTMB)

source("../data/recdata.R")

objfun <- function(par)
{
  Rmax <- exp(par[["logRmax"]])
  S50 <- exp(par[["logS50"]])
  sigma <- exp(par[["logSigma"]])

  R <- recdata$R
  S <- recdata$S

  Rhat <- Rmax * S / (S + S50)
  nll <- -sum(dnorm(log(R), log(Rhat), sigma, TRUE))
  ADREPORT(log(Rhat))
  nll
}

par <- list(logRmax=0, logS50=0, logSigma=0)
model <- MakeADFun(objfun, par, silent=TRUE)
fit <- nlminb(model$par, model$fn, model$gr)
rep <- summary(sdreport(model))

hessian <- model$he()
SE <- sqrt(diag(solve(hessian)))
cor <- round(solve(hessian) / (SE %o% SE), 3)

rep
cor
fit
