## ----setup, include=FALSE-----------------------------------------------------
library(nimble)
has_compareMCMCs <- require(compareMCMCs)
if(!has_compareMCMCs)
  warning("This Rmd file uses package compareMCMCs from github.  Sections using it will be skipped.")
has_rjags <- require(rjags)
if(!has_rjags)
  warning("This Rmd file uses package rjags.  Sections using it will be skipped.")
doComparisons <- TRUE
runExercise <- FALSE


## ----load_Ecervi_example------------------------------------------------------
source(file.path("..", "examples", "DeerEcervi", "load_DeerEcervi.R"), 
       chdir = TRUE)


## ---- fig.cap='', fig.width=12, fig.height=5----------------------------------
set.seed(123)
model <- nimbleModel(DEcode, constants = DEconstants, inits = DEinits_vals, data = DEdata)
cmodel <- compileNimble(model)
mcmcConf <- configureMCMC(model)
mcmcConf$addMonitors('farm_effect')
mcmc <- buildMCMC(mcmcConf)
cmcmc <- compileNimble(mcmc, project = model)
system.time(DEsamples <- runMCMC(cmcmc, niter = 5000))

# change to c(3,5) in live demo
par(mfrow = c(1,5))
ts.plot(DEsamples[ , 'farm_sd'], main = 'farm sd')
ts.plot(DEsamples[ , 'farm_effect[1]'], main = 'farm effect 1')
ts.plot(DEsamples[ , 'sex_int[1]'], main = 'sex intercept 1')
ts.plot(DEsamples[ , 'sex_int[2]'], main = 'sex intercept 2')
ts.plot(DEsamples[ , 'length_coef[1]'], main = 'length slope 1')


## -----------------------------------------------------------------------------
DEcodeCtrRE <- nimbleCode({
  sex_int[2] ~ dnorm(0, sd = 1000)     # offset for sex 2
  sex_int[1] <- 0                      # replaced by 'mu'
  for(i in 1:2) {
    # Priors for ntercepts and length coefficients for sex = 1,2
    length_coef[i] ~ dnorm(0, sd = 1000)
  }
  
  # Priors for farm random effects and their standard deviation.
  farm_sd ~ dunif(0, 20)
  for(i in 1:num_farms) {
    farm_effect[i] ~ dnorm(mu, sd = farm_sd)  # centered random effects
  }
  mu ~ dnorm(0, sd = 100)
  
  # logit link and Bernoulli data probabilities
  for(i in 1:num_animals) {
    logit(disease_probability[i]) <- 
      sex_int[ sex[i] ] +
      length_coef[ sex[i] ]*length[i] +
      farm_effect[ farm_ids[i] ]
    Ecervi_01[i] ~ dbern(disease_probability[i])
  }
})

set.seed(123)
modelCtrRE = nimbleModel(DEcodeCtrRE, constants = DEconstants, inits = c(DEinits_vals, mu = 0), data = DEdata)


## ---- fig.cap='', fig.width=12, fig.height=5----------------------------------
cmodelCtrRE <- compileNimble(modelCtrRE)
mcmcConfCtrRE <- configureMCMC(modelCtrRE)
mcmcConfCtrRE$addMonitors('farm_effect')
mcmcCtrRE <- buildMCMC(mcmcConfCtrRE)
cmcmcCtrRE <- compileNimble(mcmcCtrRE, project = modelCtrRE)
system.time(DEsamplesCtrRE <- runMCMC(cmcmcCtrRE, niter = 5000))

par(mfrow = c(1,5))
ts.plot(DEsamplesCtrRE[ , 'farm_sd'], main = 'farm sd')
ts.plot(DEsamplesCtrRE[ , 'farm_effect[1]'], main = 'farm effect 1')
ts.plot(DEsamplesCtrRE[ , 'mu'], main = 'intercept (random effect mean)')
ts.plot(DEsamplesCtrRE[ , 'sex_int[2]'], main = 'sex 2 offset')
ts.plot(DEsamplesCtrRE[ , 'length_coef[1]'], main = 'length slope 1')


## -----------------------------------------------------------------------------
DEconstants_uncLen <- DEconstants # Artificially un-center Length
DEconstants_uncLen$length <- DeerEcervi$Length 

set.seed(123)
DEinits_vals_uncLen <- DEinits_vals
DEinits_vals_uncLen$sex_int <- c(-8, -8)

modelUncLen = nimbleModel(DEcode, constants = DEconstants_uncLen, inits = DEinits_vals_uncLen, data = DEdata)


## ---- fig.cap='', fig.width=12, fig.height=5----------------------------------
cmodelUncLen <- compileNimble(modelUncLen)
mcmcConfUncLen <- configureMCMC(modelUncLen)
mcmcConfUncLen$addMonitors('farm_effect')
mcmcUncLen <- buildMCMC(mcmcConfUncLen)
cmcmcUncLen <- compileNimble(mcmcUncLen, project = modelUncLen)
system.time(DEsamplesUncLen <- runMCMC(cmcmcUncLen, niter = 5000))

par(mfrow = c(1,5))
ts.plot(DEsamplesUncLen[ , 'farm_sd'], main = 'farm sd')
ts.plot(DEsamplesUncLen[ , 'farm_effect[1]'], main = 'farm effect 1')
ts.plot(DEsamplesUncLen[ , 'sex_int[1]'], main = 'sex intercept 1')
ts.plot(DEsamplesUncLen[ , 'sex_int[2]'], main = 'sex intercept 2')
ts.plot(DEsamplesUncLen[ , 'length_coef[1]'], main = 'length slope 1')


## ---- fig.cap='', fig.width=12, fig.height=5----------------------------------
modelBlock = nimbleModel(DEcode, constants = DEconstants_uncLen, inits = DEinits_vals_uncLen, data = DEdata)
cmodelBlock <- compileNimble(modelBlock)
mcmcConfBlock <- configureMCMC(modelBlock)
mcmcConfBlock$removeSamplers(c('sex_int','length_coef'))

# Add RW_block samplers, modifying adaptation behavior.
mcmcConfBlock$addSampler(target = c('sex_int[1]', 'length_coef[1]'),
                 type = "RW_block",
                 control = list(propCov = diag(c(.1, .01)), adaptInterval = 20, 
                                adaptFactorExponent = 0.25))
mcmcConfBlock$addSampler(target = c('sex_int[2]', 'length_coef[2]'),
                 type = "RW_block",
                 control = list(propCov = diag(c(.1, .01)), adaptInterval = 20, 
                                adaptFactorExponent = 0.25))
mcmcConfBlock$addMonitors('farm_effect')
mcmcBlock <- buildMCMC(mcmcConfBlock)
cmcmcBlock <- compileNimble(mcmcBlock, project = modelBlock)
system.time(DEsamplesBlock <- runMCMC(cmcmcBlock, niter = 5000))

par(mfrow = c(1,5))
ts.plot(DEsamplesBlock[ , 'farm_sd'], main = 'farm sd')
ts.plot(DEsamplesBlock[ , 'farm_effect[1]'], main = 'farm effect 1')
ts.plot(DEsamplesBlock[ , 'sex_int[1]'], main = 'sex intercept 1')
ts.plot(DEsamplesBlock[ , 'sex_int[2]'], main = 'sex intercept 2')
ts.plot(DEsamplesBlock[ , 'length_coef[1]'], main = 'length slope 1')






## ---- include=FALSE-----------------------------------------------------------
mcmcResults_jags_uncLen <- list()




## ---- echo=FALSE, val=doComparisons-------------------------------------------
mcmcResults_ctrRE <- c(mcmcResults_nimble_ctrRE)
mcmcResults_uncLen <- c(mcmcResults_nimble_uncLen,
                 mcmcResults_jags_uncLen) ## These are lists of MCMCresult objects

make_MCMC_comparison_pages(mcmcResults_nimble_ctrRE, 
                           modelName = "orig_deer_ecervi_ctrRE_results")
make_MCMC_comparison_pages(mcmcResults_uncLen, modelName = "orig_deer_ecervi_uncLen_results")


## ---- eval=FALSE--------------------------------------------------------------
## # Run this code to generate your own results
## mcmcResults_ctrRE <- c(mcmcResults_nimble_ctrRE)
## mcmcResults_uncLen <- c(mcmcResults_nimble_uncLen,
##                  mcmcResults_jags_uncLen) ## These are lists of MCMCresult objects
## 
## make_MCMC_comparison_pages(mcmcResults_nimble_ctrRE,
##                            modelName = "deer_ecervi_ctrRE_results")
## make_MCMC_comparison_pages(mcmcResults_uncLen, modelName = "deer_ecervi_uncLen_results")


## -----------------------------------------------------------------------------
code_heavy <- nimbleCode({
  for(t in 1:n) 
    y[t] ~ dnorm(x[t], sd = sigma)
  for(t in 2:n) {
    x[t] <- x[t-1] + eps[t-1]
    eps[t] ~ dnorm(0, sd = omega)
  }
})


## -----------------------------------------------------------------------------
code_light <- nimbleCode({
  for(t in 1:n) 
    y[t] ~ dnorm(x[t], sd = sigma)
  for(t in 2:n)
    x[t] ~ dnorm(x[t-1], sd = omega)
})


## -----------------------------------------------------------------------------
n <- 20
m_heavy <- nimbleModel(code_heavy, 
                       data = list(y = rnorm(n)), 
                       constants = list(n = n))
m_light <- nimbleModel(code_light, 
                       data = list(y = rnorm(n)), 
                       constants = list(n = n))


## -----------------------------------------------------------------------------
m_heavy$getDependencies('eps[18]')
m_light$getDependencies('x[18]')

m_heavy$getDependencies('eps[1]')
m_light$getDependencies('x[1]')


## -----------------------------------------------------------------------------
code <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 100)
  predicted.y[1:4] <- intercept + slope * x[1:4] # vectorized node
  for(i in 1:4) {
    y[i] ~ dnorm(predicted.y[i], sd = sigma)
  }
})
model <- nimbleModel(code, data = list(y = rnorm(4)))

model$getDependencies('slope')


## -----------------------------------------------------------------------------
code <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 100)
  predicted.y[1:4] <- intercept + slope * x[1:4] # vectorized node
  for(i in 1:4) {
    x[i] ~ dnorm(0, 1)   # scalar random effects
    y[i] ~ dnorm(predicted.y[i], sd = sigma)
  }
})
model <- nimbleModel(code, data = list(y = rnorm(4)))

model$getDependencies('x[2]')


## -----------------------------------------------------------------------------
DEcodeFlex <- nimbleCode({
  sex_int[1] <- 0    # constraint to allow mixture to have non-zero mean
  sex_int[2] ~ dnorm(0, sd = 1000)
  for(i in 1:2) {
    # Priors for intercepts and length coefficients for sex = 1,2
    length_coef[i] ~ dnorm(0, sd = 1000)
  }
  
  # Priors for farm random effects
  # 'Manual' inclusion of bivariate normal mixture
  for(i in 1:num_farms) {
    farm_effect[i] ~ dnorm(mu[ind[i]+1], sd = sigma[ind[i]+1])
    ind[i] ~ dbern(pi)
  }
  for(i in 1:2) {
    mu[i] ~ dnorm(0, sd = 1000)
    sigma[i] ~ dunif(0, 20)
  }
  pi ~ dbeta(1, 1)   # same as dunif(0,1) but conjugacy will be detected
  
  # logit link and Bernoulli data probabilities
  for(i in 1:num_animals) {
    logit(disease_probability[i]) <- 
      sex_int[ sex[i] ] +
      length_coef[ sex[i] ]*length[i] +
      farm_effect[ farm_ids[i] ]
    Ecervi_01[i] ~ dbern(disease_probability[i])
  }
})


## -----------------------------------------------------------------------------
dnormmix2 <- nimbleFunction(
  run = function(x = double(0), prob = double(0), 
                 mean = double(1), sd = double(1), 
                 log = logical(0, default = 0)) {
    
    returnType(double(0))
    # generally we want to calculate probability (density) on a 
    # log scale, but here that won't work.
    dens <- prob     * dnorm(x, mean[1], sd[1]) + 
            (1-prob) * dnorm(x, mean[2], sd[2])  
    if(log) 
      return(log(dens)) else return(dens)
  })


## ---- include=FALSE-----------------------------------------------------------
# only needed for Rmd compilation; not needed for regular usage.
assign('dnormmix2', dnormmix2, .GlobalEnv)
# 'r' simulation function not required but included here because of Rmd compilation issues.
rnormmix2 <- nimbleFunction(
  run = function(n = integer(0), prob = double(0), 
                 mean = double(1), sd = double(1)) {
  # warning: dummy code    
  returnType(double(0))
  return(0)
})

assign('rnormmix2', rnormmix2, .GlobalEnv)


## -----------------------------------------------------------------------------
DEcodeFlexMarg <- nimbleCode({
  # Priors for intercepts and length coefficients for sex = 1,2
  sex_int[1] <- 0    # constraint to allow mixture to have non-zero mean
  sex_int[2] ~ dnorm(0, sd = 1000)
  for(i in 1:2) {
    length_coef[i] ~ dnorm(0, sd = 1000)
  }
  
  # Priors for farm random effects (centered on the 'baseline' sex)
  for(i in 1:num_farms) {
    farm_effect[i] ~ dnormmix2(pi, mu[1:2], sigma[1:2])
  }
  for(i in 1:2) {
    mu[i] ~ dnorm(0, sd = 1000)
    sigma[i] ~ dunif(0, 20)
  }
  pi ~ dbeta(1, 1)
  
  # logit link and Bernoulli data probabilities
  for(i in 1:num_animals) {
    logit(disease_probability[i]) <- 
      sex_int[ sex[i] ] +
      length_coef[ sex[i] ]*length[i] +
      farm_effect[ farm_ids[i] ]
    Ecervi_01[i] ~ dbern(disease_probability[i])
  }
})


## -----------------------------------------------------------------------------
set.seed(1)
modelFlexMarg <- nimbleModel(DEcodeFlexMarg, data = DEdata, 
                     constants = DEconstants, 
                     inits = c(DEinits_vals, list(pi = runif(1), 
                               mu = rnorm(2), sigma = rep(1, 2))))

modelFlexMarg$calculate('farm_effect')
cModelFlexMarg <- compileNimble(modelFlexMarg)
cModelFlexMarg$calculate('farm_effect')

