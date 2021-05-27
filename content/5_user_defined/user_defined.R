## ----setup, include=FALSE-----------------------------------------------------
library(nimble)


## -----------------------------------------------------------------------------
dnormmix2 <- nimbleFunction(
  run = function(x = double(), prob = double(), 
                 mean = double(1), sd = double(1), 
                 log = logical(0, default = 0)) {
    returnType(double(0))
    # generally calculate log density for numerical stability but not (directly) feasible here
    dens <- prob*dnorm(x, mean[1], sd[1]) + (1-prob)*dnorm(x, mean[2], sd[2])  
    # Note that use of distributions in nimbleFunctions uses default R parameterizations 
    # (e.g., sd here)
    if(log) 
      return(log(dens)) else return(dens)
  })


## -----------------------------------------------------------------------------
rnormmix2 <- nimbleFunction(
  run = function(n = integer(), 
                 prob = double(), mean = double(1), sd = double(1)) {
    returnType(double(0))
    
    # Use of distributions in nimbleFunctions follows R parameterizations
    ind <- rbinom(1, 1, prob) + 1  # dbern not available for compilation
    value <- rnorm(1, mean[ind], sd[ind])
    return(value)
  })


## -----------------------------------------------------------------------------
dnormmix2(x = 1.2, prob = 0.2, mean = c(0, 2), sd = c(0.3, 0.7), log = TRUE)
# We could use standard R debugging tools such as debug(dnormmix2) or inserting "browser()" into it.


## -----------------------------------------------------------------------------
c_dnormmix2 <- compileNimble(dnormmix2)
c_dnormmix2(x = 1.2, prob = 0.2, mean = c(0, 2), sd = c(0.3, 0.7), log = TRUE)


## -----------------------------------------------------------------------------
expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma  # calculate once
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho) # vectorized alternative is given later
    return(result)
  })


## ---- include=FALSE-----------------------------------------------------------
# only needed for Rmd compilation; not needed for regular usage.
assign('expcov', expcov, .GlobalEnv)


## -----------------------------------------------------------------------------
code <- nimbleCode({
  mu[1:N] <- mu0 * ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma)
  x[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  # other parts of model omitted
})


## -----------------------------------------------------------------------------
code <- nimbleCode({
  # (hyper)parameter priors
  mu0 ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)  # prior for variance components based on Gelman (2006)
  rho ~ dunif(0, 5)      # there might be a better non-informative prior for this

  # MVN normal (Gaussian process) prior
  mu[1:N] <- mu0 * ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma)
  x[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  
  # likelihood for count data (e.g., disease mapping)
  for(i in 1:N) {
    lambda[i] <- expected[i] * exp(x[i])
    y[i] ~ dpois(lambda[i])
  }
})

N <- 134
dists <- as.matrix(dist(runif(N)))
model <- nimbleModel(code, constants = list(N = N, dists = dists, ones = rep(1, N)), 
                     inits = list(rho = 1, sigma = 1, mu0 = 0))
deps <- model$getDependencies(c('rho','mu0','sigma'), self = FALSE)
deps
model$simulate(deps)  # may be a bit slow uncompiled given the nested looping
range(model$x)


## -----------------------------------------------------------------------------
# Verify node name
covNode <- model$expandNodeNames("cov")
covNode


## ---- eval=FALSE--------------------------------------------------------------
## debug(expcov)


## -----------------------------------------------------------------------------
model$calculate("cov[1:134, 1:134]")


## -----------------------------------------------------------------------------
expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    result <- sigma*sigma * exp(-dists / rho)
    return(result)
  })


## -----------------------------------------------------------------------------
code <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 100)
  predicted.y[1:4] <- intercept + slope * x[1:4] # vectorized node
  y[1:4] ~ dnorm_vec(predicted.y[1:4], sd = sigma)
})


## ---- eval=FALSE--------------------------------------------------------------
## my_dnorm <- nimbleFunction(
##   run = function(x = double(0),
##                  mean = double(0), sd = double(0),
##                  log = logical(0, default = 0)) {
##     returnType(double(0))
##     logdens <- dnorm(x, mean, sd, log = TRUE)
##     if(log) return(logdens) else return(exp(dens))
##   })

