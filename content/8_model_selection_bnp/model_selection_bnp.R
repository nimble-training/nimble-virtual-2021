## ----chunksetup, include=FALSE------------------------------------------------
library(nimble)


## ---- avandia-view------------------------------------------------------------
dat <- read.csv('../examples/avandia.csv')
head(dat)


## ---- avandia-setup-----------------------------------------------------------
dat <- dat[-49, ]   # This study is rather different than the others.

x <- dat$controlMI
n <- dat$nControl
y <- dat$avandiaMI
m <- dat$nAvandia

nStudies <- nrow(dat)
data <- list(x = x, y = y)
constants = list(n = n, m = m, nStudies = nStudies)
inits = list(theta = 0, mu = 0, tau = 1, gamma = rnorm(nStudies))

codeParam <- nimbleCode({
  for(i in 1:nStudies) {
    y[i] ~ dbin(size = m[i], prob = q[i]) # avandia MIs
    x[i] ~ dbin(size = n[i], prob = p[i]) # control MIs
    q[i] <- expit(theta + gamma[i])       # Avandia log-odds
    p[i] <- expit(gamma[i])               # control log-odds
    gamma[i] ~ dnorm(mu, sd = tau)        # study effects
  }
  theta ~ dflat()        # effect of Avandia
  # random effects hyperparameters
  mu ~ dflat()
  tau ~ dunif(0, 100)
})


## ---- mcmc, fig.cap='', fig.width=12, fig.height=5----------------------------
samples <- nimbleMCMC(code = codeParam, data = data, inits = inits, 
                      constants = constants, monitors = c("mu", "tau", "theta", "gamma"),
                      thin = 10, niter = 11000, nburnin = 1000, nchains = 1, setSeed = TRUE)
gammaCols <- grep('gamma', colnames(samples))

par(mfrow = c(1, 4))
ts.plot(samples[ , 'theta'], xlab = 'iteration', ylab = expression(theta))
hist(samples[ , 'theta'], xlab = expression(theta), main = 'effect of Avandia')
gammaMn <- colMeans(samples[ , gammaCols])
hist(gammaMn, xlab = 'posterior means of random effects', main = 'random effects distribution')
hist(samples[500, gammaCols], xlab = 'single draw of random effects',
                   main = 'random effects distribution')


## ---- mcmc-output-------------------------------------------------------------
# Posterior probability of a positive coefficient
mean(samples[ , 'theta'] > 0)
# What about a substantively significant effect
cutoff <- 0.15    # Suppose theta > 0.15 is deemed non-negligible
mean(samples[ , 'theta'] > cutoff)


## ---- waic--------------------------------------------------------------------
outputFull <- nimbleMCMC(code = codeParam, data = data, inits = inits,
                      constants = constants, monitors = c("mu", "tau", "gamma", "theta"), 
                      thin = 10, niter = 11000, nburnin = 1000, nchains = 1, 
                      setSeed = TRUE, WAIC = TRUE)

codeParamReduced <- nimbleCode({
    for(i in 1:nStudies) {
        y[i] ~ dbin(size = m[i], prob = q[i]) # avandia MIs
        x[i] ~ dbin(size = n[i], prob = p[i]) # control MIs
        q[i] <- expit(gamma[i])               # Avandia arms log-odds; no Avandia effect
        p[i] <- expit(gamma[i])               # control log-odds
        gamma[i] ~ dnorm(mu, sd = tau)        # study effects
    }
    # random effects hyperparameters
    mu ~ dflat()
    tau ~ dunif(0, 100)
})

outputReduced <- nimbleMCMC(code = codeParamReduced, data = data, inits = inits,
                      constants = constants, monitors = c("mu", "tau", "gamma"),
                      thin = 10, niter = 11000, nburnin = 1000, nchains = 1, 
                      setSeed = TRUE, WAIC = TRUE)


## -----------------------------------------------------------------------------
outputFull$WAIC
outputReduced$WAIC


## ---- meta-rjmcmc-setup-------------------------------------------------------
codeParam2 <- nimbleCode({
  for(i in 1:nObs) {
    full_y[i] ~ dbin(size = full_n[i], prob = p[i])
    p[i] <- expit(theta*avandia[i] + gamma[study[i]])       #log-odds
  }
  for(i in 1:nStudies)
    gamma[i] ~ dnorm(mu, sd = tau)        # study effects
  theta ~ dflat()        # effect of Avandia
  # random effects hyperparameters
  mu ~ dflat()
  tau ~ dunif(0, 100)
})

full_y <- c(dat$controlMI, dat$avandiaMI)
avandia <- c(rep(0, nrow(dat)), rep(1, nrow(dat)))
full_n <-  c(dat$nControl, dat$nAvandia)

nObs <- 2*nrow(dat)
data2 <- list(full_y = full_y)
constants2 = list(full_n = full_n, study = rep(1:nStudies, 2), 
                 nObs = nObs, avandia = avandia, nStudies = nStudies)

model2 <- nimbleModel(code = codeParam2, data = data2, inits = inits, constants = constants2)
cModel2 <- compileNimble(model2)


## -----------------------------------------------------------------------------
conf2 <- configureMCMC(model2, monitors = c("mu", "tau", "theta", "gamma"))
configureRJ(conf2,
            targetNodes = 'theta',
            priorProb = 0.5,
            control = list(mean = 0, scale = 1))
mcmc2 <- buildMCMC(conf2)

cmcmc2 <- compileNimble(mcmc2, project = model2)
resultsVarSel <- runMCMC(cmcmc2, niter = 11000, nburnin = 1000, thin = 10, setSeed = 1)


## ---- results, fig.width=8, fig.height=5, fig.cap=''--------------------------
par(mfrow = c(1,2))
ts.plot(resultsVarSel[ , 'theta'] != 0, xlab = 'iteration', ylab = 'theta presence',
                    main = 'Avandia effect presence')
ts.plot(resultsVarSel[ , 'theta'], xlab = 'iterations', ylab = 'theta',
               main = 'Avandia effect')

## posterior probability of inclusion    
mean(resultsVarSel[ , 'theta'] != 0)  


## ---- meta-bnp----------------------------------------------------------------
codeBNP <- nimbleCode({
  for(i in 1:nStudies) {
    y[i] ~ dbin(size = m[i], prob = q[i]) # avandia MIs
    x[i] ~ dbin(size = n[i], prob = p[i]) # control MIs
    q[i] <- expit(theta + gamma[i])       # Avandia log-odds
    p[i] <- expit(gamma[i])               # control log-odds
    
    # Dirichlet process prior for random effects
    gamma[i] ~ dnorm(mu[i], var = tau[i]) # random effects (from mixture)
    mu[i] <- muTilde[xi[i]]               # mean for component assigned to i'th study
    tau[i] <- tauTilde[xi[i]]             # variance for component assigned to i'th study
  }
  # mixture component parameters drawn from base measures
  # should think carefully about priors here
  for(i in 1:nStudies) {
    muTilde[i] ~ dnorm(-6, sd = 1)  # based on parametric fit (slightly cheating...)
    tauTilde[i] ~ dinvgamma(2, 1)
  }
  # CRP for clustering studies to mixture components
  xi[1:nStudies] ~ dCRP(conc, size = nStudies)
  # hyperparameters
  conc ~ dgamma(1, 1)      # 'alpha' in the CRP discussion
  theta ~ dflat()          # effect of Avandia
})


## ---- DP-MCMC, fig.cap='', fig.width=12, fig.height=5-------------------------
set.seed(1)
inits <- list(gamma = rnorm(nStudies), xi = sample(1:2, nStudies, replace = TRUE),
              conc = 1, theta = 0,
              muTilde = rnorm(nStudies), tauTilde = rep(1, nStudies))

samplesBNP <- nimbleMCMC(code = codeBNP, data = data, inits = inits,
                         constants = constants,
                         monitors = c("theta", "gamma", "conc", "xi"),
                         thin = 10, niter = 11000, nburnin = 1000, nchains = 1)

gammaCols <- grep('gamma', colnames(samplesBNP))
xiCols <- grep('xi', colnames(samplesBNP))

par(mfrow = c(1,5))
ts.plot(samplesBNP[ , 'theta'], xlab = 'iteration', ylab = expression(theta))
hist(samplesBNP[ , 'theta'], xlab = expression(theta), main = 'effect of Avandia')
gammaMn <- colMeans(samplesBNP[ , gammaCols])
hist(gammaMn, xlab = 'posterior means of random effects',
     main = 'random effects distribution')
hist(samplesBNP[1000, gammaCols], xlab = 'single draw of random effects',
     main = 'random effects distribution')

# How many mixture components are inferred?
xiRes <- samplesBNP[ , xiCols]
nGrps <- apply(xiRes, 1, function(x) length(unique(x)))
ts.plot(nGrps, xlab = 'iteration', ylab = 'number of components')


## ---- DP-samplers-------------------------------------------------------------
model <- nimbleModel(codeBNP, constants = constants, data = data, inits = inits)
conf = configureMCMC(model, print = TRUE)

