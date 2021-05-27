## ----chunksetup, include=FALSE------------------------------------------------
library(nimble)


## ---- fig.cap = '', fig.height=4, fig.width=4, echo = FALSE-------------------
locs <- cbind(c(-1,1,-1,1),c(-1,-1,1,1)) 
plot(locs, type = 'n', xlab = '', ylab = '', xaxt = 'n', yaxt = 'n')
abline(v=0); abline(v=2); abline(v=-2)
abline(h=0); abline(h=2); abline(h=-2)
text(locs[,1]/2, locs[,2]/2, as.character(1:4), cex = 3)


## ---- eval=FALSE--------------------------------------------------------------
## num <- c(2, 2, 2, 2)  # each location has two neighbors
## adj <- c(2, 3,        # neighbors of loc'n 1
## 	      1, 4,         # neighbors of loc'n 2
## 	      1, 4,         # etc.
## 	      2, 3)


## ---- fig.cap='', fig.width=12, fig.height=7----------------------------------
library(CARBayesdata, quietly = TRUE)
library(sp, quietly = TRUE)
library(spdep, quietly = TRUE)
library(classInt)
data(GGHB.IG)
data(respiratorydata)

respiratorydata_sp <- merge(x = GGHB.IG, y = respiratorydata, by = "IG", all.x = FALSE)
respiratorydata_sp <- spTransform(respiratorydata_sp,
                                  CRS("+proj=longlat +datum=WGS84 +no_defs"))

if(FALSE) { 
  # This will produce an image on top of an OpenStreetMap webpage}
  library(leaflet)
  colors <- colorNumeric(palette = "YlOrRd", domain = respiratorydata_sp@data$SMR)
  map2 <- leaflet(data=respiratorydata_sp) %>%
    addPolygons(fillColor = ~colors(SMR), color="", weight=1,
                fillOpacity = 0.7) %>%
    addLegend(pal = colours, values = respiratorydata_sp@data$SMR, opacity = 1,
              title="SMR") %>%
    addScaleBar(position="bottomleft")
  map2
}

# Instead, make a simple map in R
temp.colors<-function(n = 25) {
  m <- floor(n/2)
  blues <- hsv(h=.65, s=seq(1,0,length=m+1)[1:m])
  reds <- hsv(h=0, s=seq(1,0,length=m+1)[1:m])
  c(blues,if(n%%2!=0) "#FFFFFF", reds[m:1])
}

q <- classIntervals(respiratorydata_sp@data$SMR ,style = "fixed",
                    fixedBreaks=seq(0.3, 1.7, by = 0.1))
pal <- temp.colors(12)
col <- findColours(q, pal)
plot(respiratorydata_sp, col = col, border = 'grey', lwd = 0.4)


## -----------------------------------------------------------------------------
W_nb <- poly2nb(respiratorydata_sp, row.names =  rownames(respiratorydata_sp@data))
## Determine neighborhood/adjacency information needed for neighborhood-based CAR model
nbInfo <- nb2WB(W_nb)

# A vector of indices indicating which regions are neighbors of which.
head(nbInfo$adj, n = 30)
# A vector of weights. In this case, all weights are 1.
head(nbInfo$weights)
# A vector of length N indicating how many neighbors each region has.
# This helps map the adj vector to each region.
nbInfo$num


## -----------------------------------------------------------------------------
code <- nimbleCode({
  # priors
  beta ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)   # prior for variance components based on Gelman (2006)
  tau <- 1 / sigma^2
  # latent process
  x[1:N] ~ dcar_normal(adj[1:L], weights[1:L], num[1:N], tau, zero_mean = 0)
  # likelihood
  for(i in 1:N) {
    lambda[i] <- expected[i] * exp(beta*z[i] + x[i])
    y[i] ~ dpois(lambda[i])
  }
})

z <- respiratorydata_sp$incomedep
z <- z - mean(z)  # center for improved MCMC performance

set.seed(1)

nregions <- nrow(respiratorydata_sp)
expected <- respiratorydata_sp$expected
y <- respiratorydata_sp$observed
constants <- list(N = nregions, L = length(nbInfo$adj), 
                  adj = nbInfo$adj, weights = nbInfo$weights, num = nbInfo$num,
                  z = z, expected = expected)
data <- list(y = y)
inits <- list(beta = 0, sigma = 1, x = rnorm(nregions))

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)

conf <- configureMCMC(model, monitors = c('beta', 'sigma', 'x'))
conf$printSamplers()

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)
samples <- runMCMC(cMCMC, niter = 5000, nburnin = 1000, thin = 10)


## ---- fig.cap='', fig.height=5, fig.width=12----------------------------------
par(mfrow = c(1,4))
ts.plot(samples[ , 'sigma'], main = 'sigma')
ts.plot(samples[ , 'x[5]'], main = 'x[5]')
ts.plot(samples[ , 'x[50]'], main = 'x[50]')
ts.plot(samples[ , 'beta'], main = 'beta')


## ---- fig.caption='', fig.width=12, fig.height=7------------------------------
xCols <- grep('^x\\[', colnames(samples))
xEstCAR <- colMeans(samples[ , xCols])

q <- classIntervals(xEstCAR,style = "fixed",fixedBreaks=seq(-0.8, 0.8, by = 0.1))
pal <- temp.colors(16)
col <- findColours(q, pal)

par(mfrow = c(1,1))
plot(respiratorydata_sp, col = col, border = 'grey', lwd = 0.4)


## -----------------------------------------------------------------------------
expcov <- nimbleFunction(     
  run = function(dists = double(2), rho = double(0), sigma = double(0)) {
    returnType(double(2))
    n <- dim(dists)[1]
    result <- matrix(nrow = n, ncol = n, init = FALSE)
    sigma2 <- sigma*sigma   # calculate once
    for(i in 1:n)
      for(j in 1:n)
        result[i, j] <- sigma2*exp(-dists[i,j]/rho)
    return(result)
  })


## ---- include=FALSE-----------------------------------------------------------
# only needed for Rmd compilation; not needed for regular usage.
assign('expcov', expcov, .GlobalEnv)


## -----------------------------------------------------------------------------
code <- nimbleCode({
  mu0 ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)  # prior for variance components based on Gelman (2006)
  rho ~ dunif(0, 5)
  beta ~ dnorm(0, sd = 100)
  
  # latent spatial process
  mu[1:N] <- mu0*ones[1:N]
  cov[1:N, 1:N] <- expcov(dists[1:N, 1:N], rho, sigma)
  x[1:N] ~ dmnorm(mu[1:N], cov = cov[1:N, 1:N])
  # likelihood
  for(i in 1:N) {
    lambda[i] <- expected[i] * exp(beta*z[i] + x[i])
    y[i] ~ dpois(lambda[i])
  }
})

locs <-  as.matrix(respiratorydata_sp@data[ , c('easting', 'northing')]) / 1e5
dists <- as.matrix(dist(locs))
dists <- dists / max(dists)  # normalize to max distance of 1

constants <- list(N = nregions, dists = dists, ones = rep(1, nregions),
                            z = z, expected = expected)
data <- list(y = y)
inits <- list(beta = 0, mu0 = 0, sigma = 1, rho = 0.2)

set.seed(1)

## setup initial spatially-correlated latent process values
inits$cov <- expcov(dists, inits$rho, inits$sigma)
inits$x <-  t(chol(inits$cov)) %*% rnorm(nregions)
inits$x <- inits$x[ , 1]  # so can give nimble a vector rather than one-column matrix

model <- nimbleModel(code, constants = constants, data = data, inits = inits)
cModel <- compileNimble(model)


## -----------------------------------------------------------------------------
conf <- configureMCMC(model)
conf$addMonitors('x')
conf$removeSamplers('x[1:134]')
## Changing a tunable parameter in the adaptation of RW_block makes a big difference.
conf$addSampler('x[1:134]', 'RW_block', control = list(adaptFactorExponent = 0.25))

MCMC <- buildMCMC(conf)
cMCMC <- compileNimble(MCMC, project = cModel)
samples <- runMCMC(cMCMC, niter = 10000, nburnin = 5000, thin = 5)


## ---- fig.cap='', fig.width=10, fig.height=6----------------------------------
par(mfrow = c(2,3))
ts.plot(samples[ , 'rho'], main = 'rho')
ts.plot(samples[ , 'sigma'], main = 'sigma')
ts.plot(samples[ , 'beta'], main = 'beta')
ts.plot(samples[ , 'x[5]'], main = 'x[5]')
ts.plot(samples[ , 'x[50]'], main = 'x[50]')
ts.plot(samples[ , 'x[100]'], main = 'x[100]')


## ---- fig.cap='', fig.width=12, fig.height=7----------------------------------
xCols <- grep('^x\\[', colnames(samples))
xEstGP <- colMeans(samples[ , xCols])

q <- classIntervals(xEstGP, style = "fixed", fixedBreaks=seq(-1.2, 1.2, by = 0.2))
pal <- temp.colors(12)
col <- findColours(q, pal)

par(mfrow = c(1,1))
plot(respiratorydata_sp, col = col, border = 'grey', lwd = 0.4)


## -----------------------------------------------------------------------------
K <- 40 # number of spline basis functions
# fit a simple GAM to get a sense for number of basis functions; K = 40 might be too small...
if(FALSE) {
  mod <- gam(y ~ log(expected) + z + s(locs, k = K), family = 'poisson')
  summary(mod)
  fit <- predict(mod, type = 'terms')[,3]
}

# ignore z to get basis info, but that may affect initial values.
out <- mgcv::jagam(y ~ log(expected) + s(locs, k = K) -1, family = 'poisson', 
                   na.action = na.pass, file = 'blank.jags')
B <- out$jags.data$X[ , -1]
K <- ncol(B)   # actually only 39 basis functions
## two components of precision matrix to ensure propriety
## of prior on coeffs (see Wood (2016) J Statistical Software)
Q1 <- out$jags.data$S1[ , 1:K]
Q2 <- out$jags.data$S1[ , (K+1):(2*K)]  


## -----------------------------------------------------------------------------
codeS <- nimbleCode({
  mu0 ~ dflat()
  beta_z ~ dnorm(0, sd = 100)
  for(i in 1:2) 
    lambda[i]  ~ dgamma(.05, .005)  # based on jagam defaults
  
  # latent process is product of (known) spline basis and (unknown) coefficients
  # two components to ensure propriety of prior on coeffs
  prec[1:K, 1:K] <- lambda[1] * Q1[1:K, 1:K] + lambda[2] * Q2[1:K, 1:K]
  beta_x[1:K] ~ dmnorm(zeros[1:K], prec[1:K, 1:K])
  x[1:N] <-  (B[1:N, 1:K] %*% beta_x[1:K])[,1]
  
  # likelihood
  for(i in 1:N) {
    mu[i] <- expected[i] * exp(mu0 + beta_z*z[i] + x[i])
    y[i] ~ dpois(mu[i])
  }
})

constantsS <- list(Q1 = Q1, Q2 = Q2, B = B, K = K, N = nregions, z = z, zeros = rep(0, K), 
                            expected = expected)
dataS <- list(y = y)
initsS <- list(beta_z = 0)
## Use inits from jagam, but note that jagam call did not use 'z' or have an intercept
initsS$beta_x <-  out$jags.ini$b[2:(K+1)]
initsS$mu0 <- 0
initsS$lambda <- out$jags.ini$lambda

set.seed(1)
modelS <- nimbleModel(codeS, constants = constantsS, data = dataS, inits = initsS)
cModelS <- compileNimble(modelS)


## -----------------------------------------------------------------------------
confS <- configureMCMC(modelS)
confS$removeSamplers('beta_x[1:39]')
# Using non-default value of 0.25 for a tuning parameter of the adaptive block sampler 
# helps a lot compared to default NIMBLE or to JAGS
confS$addSampler('beta_x[1:39]', 'RW_block', control = list(adaptFactorExponent = 0.25))
confS$addMonitors('x', 'beta_x')

mcmcS <- buildMCMC(confS)
cmcmcS <- compileNimble(mcmcS, project = cModelS)
samplesS <- runMCMC(cmcmcS, niter = 100000, nburnin = 20000, thin = 80)


## ---- fig.cap='', fig.width=12, fig.height=6----------------------------------
par(mfrow = c(2,3))
ts.plot(samplesS[ , 'lambda[1]'], main = 'lambda[1]')
ts.plot(samplesS[ , 'lambda[2]'], main = 'lambda[2]')
ts.plot(samplesS[ , 'beta_z'], main = 'beta_z')
ts.plot(samplesS[ , 'x[5]'], main = 'x[5]')
ts.plot(samplesS[ , 'x[50]'], main = 'x[50]')
ts.plot(samplesS[ , 'x[100]'], main = 'x[100]')


## ---- fig.cap='', fig.width=12, fig.height=7----------------------------------
xCols <- grep('^x\\[', colnames(samplesS))
xEstSpl <- colMeans(samplesS[ , xCols])

q <- classIntervals(xEstSpl, style = "fixed", fixedBreaks=seq(-0.2, 0.2, by = 0.04))
pal <- temp.colors(10)
col <- findColours(q, pal)
par(mfrow = c(1,1))
plot(respiratorydata_sp, col = col, border = 'grey', lwd = 0.4)


## ---- eval=FALSE--------------------------------------------------------------
## newlocs <- rbind(c(2.6, 6.7), c(2.61, 6.69), c(2.59, 6.69))
## dist11 <- fields::rdist(locs)
## dist21 <- fields::rdist(newlocs, locs)
## dist22 <- fields::rdist(newlocs)


## ---- eval=FALSE--------------------------------------------------------------
## sample_xstar <- nimbleFunction({
##   # need types for inputs
##   run = function(x = double(1), mu = double(1), muNew = double(1), sigma = double(0), rho = double(0),
##                  dist11 = double(2), dist22 = double(2), dist21 = double(2)) {
##     returnType(double(1))
##     n <- length(muNew)
##     sigma2 <- sigma*sigma
##     C22 <- sigma2 * exp(-dist22 / rho)
##     C11 <- sigma2 * exp(-dist11 / rho)
##     C21 <- sigma2 * exp(-dist21 / rho)
##     # Note that this could be made a bit more efficient by using the Cholesky
##     # decomposition rather than solve().
##     xstar <- munew + (C21 %*% solve(C11, x - mu))[,1]
##     xstar <- xstar + (t(chol(C22 - C21 %*% solve(C11, t(C21)))) %*% rnorm(n))[,1]
##     return(xstar)
##   }
## })
## 
## get_samples <- nimbleFunction(
##   run = function(samples = double(2), z = double(1), zNew = double(1),
##                  dist11 = double(2), dist22 = double(2), dist21 = double(2)) {
##     returnType(double(2))
##     m <- dim(samples)[1]
##     nstar <- dim(dst21)[2]
##     output <- matrix(0, nrow = m, ncol = nstar)
##     for(i in 1:m) {
##       # Extract parameter values from the input matrix based on numeric column indexes.
##       # Inelegant because hard-coded based on looking at column names of samples matrix.
##       mu0 <- output[i, 2]
##       beta <- output[i, 1]
##       rho <- output[i, 3]
##       sigma <- output[i, 4]
##       x <- output[i, 5:138]
##       mu <- mu0 + beta*z
##       muNew <- mu0 + beta*zNew
##       # Get other parameters
##       output[i, ] <- sample_xstar(x, mu, muNew, sigma, rho, dist11, dist22, dist21)
##     }
##     return(output)
##   }
## )
## 
## ## Now compare uncompiled and compiled speed of carrying out the sampling.
## xstar_samples <- get_samples(samples, z, zNew, dist11, dist22, dist21)
## cget_samples <- compileNimble(get_samples)
## xstar_samples2 <- cget_samples(samples, z, Znew, dist11, dist22, dist21)

