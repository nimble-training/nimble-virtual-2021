## ----loadLibs, include=FALSE-------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
doMCMC <- TRUE
library(nimble)
if(!require(mvtnorm))
  warning("This Rmd file needs package mvtnorm.")
if(!require(coda))
  warning("This Rmd file needs package coda")
has_compareMCMCs <- require(compareMCMCs)
if(!has_compareMCMCs)
  warning("This Rmd file uses package compareMCMCs from github.  Sections using it will be skipped.")
has_rjags <- require(rjags)
if(!has_rjags)
  warning("This Rmd file uses package rjags.  Sections using it will be skipped.")
doComparisons <- FALSE


## ----setup-------------------------------------------------------------------------------------------------------------
source(file.path("..", "examples", "DeerEcervi", "load_DeerEcervi.R"), 
       chdir = TRUE)
set.seed(123)
DEmodel <- nimbleModel(DEcode,
                        constants = DEconstants,
                        data = list(Ecervi_01 = DeerEcervi$Ecervi_01),
                        inits = DEinits())


## ----lmModel, include=FALSE--------------------------------------------------------------------------------------------
lmCode <- nimbleCode({
  sigma ~ dunif(0, 20)
  intercept ~ dnorm(0, sd = 100)
  slope ~ dnorm(0, sd = 100)
  for(i in 1:n) {
    y[i] ~ dnorm(intercept + slope * x[i], sd = sigma)
  }
})
n <- 5
lmModel <- nimbleModel(lmCode, 
                       constants = list(n = n))
lmModel$slope <- 0.6
lmModel$intercept <- 10
lmModel$sigma <- 0.2
lmModel$x <- seq(0, 1, length = n)
lmModel$calculate()
set.seed(0)
lmModel$simulate('y')
lmModel$setData('y')
# {
#   plot(lmModel$x, lmModel$y, pch = 19)
#   abline(lm(lmModel$y ~ lmModel$x))
# }


## ----show-lm-data------------------------------------------------------------------------------------------------------
data.frame(x = lmModel$x, y = lmModel$y)


## ----calc-lm-posterior, include=FALSE----------------------------------------------------------------------------------
log_posterior_numerator <- function(params) {
  lmModel$intercept <- params[1]
  lmModel$slope <- params[2]
  lmModel$calculate()
}
optim_map <- optim(c(10, 0.8), log_posterior_numerator, control = list(fnscale = -1))
optim_map$par
lmFit <- lm(lmModel$y ~ lmModel$x)
lmCoef <- coefficients(summary(lm(lmModel$y ~ lmModel$x))) ## Check that they match mle.
lmCoef
## Make a grid +/- 3 standard errors around the MLE
intercept_grid <- lmCoef['(Intercept)', 'Estimate'] +
  lmCoef['(Intercept)', 'Std. Error'] * seq(-3, 3, length = 21)
slope_grid <- lmCoef['lmModel$x', 'Estimate'] +
  lmCoef['lmModel$x', 'Std. Error'] * seq(-3, 3, length = 21)
llh_surface <- matrix(0, nrow = length(intercept_grid), 
                      ncol = length(slope_grid))
for(i in seq_along(intercept_grid))
  for(j in seq_along(slope_grid))
    llh_surface[i, j] <- log_posterior_numerator(c(intercept_grid[i], slope_grid[j]))
library(mvtnorm)
## we will "cheat" and use the mle
samples <- rmvnorm(1000, mean = lmCoef[, "Estimate"], sigma = vcov(lmFit))


## ----contour-lm-1, echo=FALSE------------------------------------------------------------------------------------------
contour(intercept_grid, slope_grid, llh_surface, 
        levels = optim_map$value - 0.01 - 0:5,
        main = "posterior density contours",
        xlab = "intercept", ylab = "slope")


## ----contour-lm-2, echo=FALSE------------------------------------------------------------------------------------------
{
  contour(intercept_grid, slope_grid, llh_surface, 
        levels = optim_map$value - 0.01 - 0:5,
        main = "posterior log density contours",
        xlab = "intercept", ylab = "slope")
  points(samples, pch = '.', col = 'red')
}


## ----mcmcConf----------------------------------------------------------------------------------------------------------
mcmcConf <- configureMCMC(DEmodel)
mcmcConf$printSamplers()








## ----lmConditional, echo=FALSE-----------------------------------------------------------------------------------------
intercept_grid_given_slope <- seq(10.1, 10.6, length = 31)
llh_surface_given_slope <- apply(matrix(intercept_grid_given_slope), 1, 
                                 function(int) log_posterior_numerator(c(int, 0.2)))
{
  plot(intercept_grid_given_slope, exp(llh_surface_given_slope), type = 'l',
     main = "Conditional posterior density (up to a constant) for slope = 0.2",
     ylab = "Conditional posterior density (up to a constant)",
     xlab = "intercept")
  lines(c(10.4, 10.4), c(0, 0.1), col = 'red')
  legend("topleft", col = "red", legend = "current value", lty = 1)
}


## ---- echo=FALSE-------------------------------------------------------------------------------------------------------
theta1 <- seq(0.5, 5, length = 200)
targetDist <- 0.1 * dnorm(theta1, 2, 0.5)
current <- 1.3
proposalDist <- dnorm(theta1, current, sd = 0.1)
proposalDisplayScale <- max(proposalDist)/max(targetDist)
proposalDist <- proposalDist / proposalDisplayScale
proposal <- 1.5
nextTargetDist <- 0.03 * dnorm(theta1, 2.4, 0.2)
{
  plot(theta1, targetDist, type = 'l', col = 'black',
       main = "Random-walk Metropolis-Hastings",
       ylab = "Target and proposal distributions (scaled)",
       xlab = expression(theta[1]))
  points(theta1, proposalDist, type = 'l', col = 'blue')
  points(theta1, nextTargetDist, type = 'l', col = 'goldenrod')
  points(current, 0.1 * dnorm(current, 2, 0.5), pch = 19, col = 'red')
  points(proposal, 0.1 * dnorm(proposal, 2, 0.5), pch = 8, col = 'red')
  lines(c(current, current), c(0, 0.1 * dnorm(current, 2, 0.5)), col = 'red')
  lines(c(proposal, proposal), c(0, 0.1 * dnorm(proposal, 2, 0.5)), col = 'red')
  legend("topright", lty = c(1,1,0,0, 1), 
         pch = c(NA, NA, 19, 8, NA), 
         col = c('black','blue','red','red', 'goldenrod'),
         legend = c('target distribution', 'proposal distribution (scaled)', 'current value', 'proposal value', 'next iteration target distribution' ))
}


## ---- echo = FALSE-----------------------------------------------------------------------------------------------------
theta1grid <- seq(0.5, 5, length = 200)
targetDist <- function(theta1) {0.1 * dnorm(theta1, 2, 0.5)}
targetDistGrid <- targetDist(theta1grid)
current <- 1.3
origCurrent <- current
update <- function(current, targetDist, mean) {
  currentF <- targetDist(current)
  u <- runif(1, 0, currentF)
  leftBound <- uniroot(function(x) (targetDist(x)-u), lower = -10, upper = mean)$root
  rightBound <- uniroot(function(x) (targetDist(x)-u), lower = mean, upper = 10)$root
  updated <- runif(1, leftBound, rightBound)
  list(lower = leftBound, upper = rightBound, u = u, currentF = currentF, updated = updated)
}
set.seed(345)
u1 <- update(current, targetDist, 2)
u2 <- update(u1$updated, targetDist, 2)
u3 <- update(u2$updated, targetDist, 2)
slicePlotter <- function(u, current) {
  points(current, targetDist(current), pch = 19, col = 'red')
  lines(c(current, current), c(0, targetDist(current)), col = 'red')
  points(current, u$u, pch = 17, col = 'purple')
  points(u$lower, u$u, pch= 3, col = 'blue')
  points(u$upper, u$u, pch = 3, col = 'blue')
  lines(c(u$lower, u$upper), c(u$u, u$u), type = 'l', col = 'blue')
  points(u$updated, u$u, pch = 8, col = 'red')
  lines(c(u$updated, u$updated), c(0, u$u), type = 'l', col = 'red')
  u$updated
}
sliceLegend <- function() {
  legend("topright",
         lty = c(1, 0, 0, 1, 0),
         legend = c('target dist', 'current value', 'vertical uniform draw', 'horizontal range', 'new value (horizontal uniform draw)'),
         pch = c(NA, 19, 17, 3, 8),
         col = c('black', 'red', 'purple','blue', 'red'))
}


## ---- echo=FALSE-------------------------------------------------------------------------------------------------------
{
  plot(theta1grid, targetDistGrid, type = 'l', col = 'black',
       main = 'slice sampler',
       ylab = "target distribution",
       xlab = expression(theta[1]))
  current <- origCurrent
  for(u in list(u1)) {
    current <- slicePlotter(u, current)
  }
  sliceLegend()
}


## ---- echo=FALSE-------------------------------------------------------------------------------------------------------
{
  plot(theta1grid, targetDistGrid, type = 'l', col = 'black',
       main = 'slice sampler',
       ylab = "target distribution",
       xlab = expression(theta[1]))
  current <- u1$updated
  for(u in list(u2)) {
    current <- slicePlotter(u, current)
  }
  sliceLegend()
}


## ---- echo=FALSE-------------------------------------------------------------------------------------------------------
{
  plot(theta1grid, targetDistGrid, type = 'l', col = 'black',
       main = 'slice sampler',
       ylab = "target distribution",
       xlab = expression(theta[1]))
  current <- u2$updated
  for(u in list(u3)) {
    current <- slicePlotter(u, current)
  }
  points(theta1grid, 0.05*dnorm(theta1grid, 3, .4), type = 'l', col = 'goldenrod')
  sliceLegend()
  legend("topleft", col = "goldenrod", pch = NA, lty = 1, legend = "new target distribution")
}




## ----------------------------------------------------------------------------------------------------------------------
nimbleMCMCdefs = list(
  nimble_RWblock = function(model) { # Input should be a model
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c('sex_int','length_coef'))
    mcmcConf$addSampler(target = c('sex_int[1]', 'length_coef[1]'),
                        type = "RW_block",
                        control = list(adaptInterval = 20, tries = 2))
    mcmcConf$addSampler(target = c('sex_int[2]', 'length_coef[2]'),
                        type = "RW_block",
                        control = list(adaptInterval = 20, tries = 2))
    mcmcConf                         # Output should be an MCMC configuration 
  },
  nimble_AFSSblock = function(model) {
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c('sex_int','length_coef'))
    mcmcConf$addSampler(target = c('sex_int[1]', 'length_coef[1]'),
                        type = "AF_slice",
                        control = list(sliceAdaptFactorInterval = 20))
    mcmcConf$addSampler(target = c('sex_int[2]', 'length_coef[2]'),
                        type = "AF_slice",
                        control = list(sliceAdaptFactorInterval = 20))
    mcmcConf                         # Output should be an MCMC configuration 
  },
  nimble_log_sigma = function(model) {
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers('farm_sd')
    mcmcConf$addSampler(target = 'farm_sd',
                        type = "RW",
                        control = list(log = TRUE))
    mcmcConf
  }
)


## ----------------------------------------------------------------------------------------------------------------------
DEcode_jags <- nimbleCode({
  for(i in 1:2) {
    length_coef[i] ~ dnorm(0, 1.0E-6) # precisions
    sex_int[i] ~ dnorm(0, 1.0E-6)
  }
  farm_sd ~ dunif(0, 20)
  farm_precision <- 1/(farm_sd*farm_sd)

  for(i in 1:num_farms) {
    farm_effect[i] ~ dnorm(0, farm_precision) # precision
  }

  for(i in 1:num_animals) {
    logit(disease_probability[i]) <-
      sex_int[ sex[i] ] +
      length_coef[ sex[i] ]*length[i] +
      farm_effect[ farm_ids[i] ]
    Ecervi_01[i] ~ dbern(disease_probability[i])
  }
})




## ---- eval=FALSE-------------------------------------------------------------------------------------------------------
## mcmcResults_nimble[['nimble_log_sigma']]$samples # output not shown
## mcmcResults_nimble[['nimble_log_sigma']]$metrics$byParameter # output not shown


## ---- include=FALSE----------------------------------------------------------------------------------------------------
mcmcResults_jags <- list()






## ---- eval=FALSE-------------------------------------------------------------------------------------------------------
## # Run this code to generate your own results
## mcmcResults <- c(mcmcResults_nimble, mcmcResults_jags) ## These are lists of MCMCresult objects
## make_MCMC_comparison_pages(mcmcResults, modelName = "deer_ecervi_mcmc_results")

