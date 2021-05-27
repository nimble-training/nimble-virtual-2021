## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
has_nimbleEcology <- require(nimbleEcology)
if(!has_nimbleEcology)
  message("This module will use nimbleEcology, which you don't have installed.")
doDHMMexample <- FALSE


## ----load_dipper--------------------------------------------------------------
dipper_example_dir <- here::here("content", "examples","Dipper")
dipper <- read.csv(file.path(dipper_example_dir,"dipper.csv"))
y <- as.matrix(dipper[ , 1:7])
first <- apply(y, 1, function(x) min(which(x !=0)))
y <- y[ first != 7, ]
y_orig <- y # Keep an "original version" for use below
head(y)


## ----dipper_basic-------------------------------------------------------------
dipper_code_basic <- nimbleCode({
  phi ~ dunif(0, 1) # survival prior
  p ~ dunif(0, 1)   # detection prior
  # likelihood
  for (i in 1:N){
    z[i,first[i]] <- 1
    for (t in (first[i]+1):T){
      z[i,t] ~ dbern(phi * z[i,t-1]) # z = 0 for dead, 1 for alive
      y[i,t] ~ dbern(p * z[i,t])     # y = 0 for not-observed, 1 for observed
    }}
  })


## ----setupInputs--------------------------------------------------------------
first <- apply(y_orig, 1, function(x) min(which(x !=0)))
dipper_constants <- list(N = nrow(y), 
                         T = ncol(y_orig), 
                         first = first)
dipper_data <- list(y = y_orig)   # 0s and 1s
zinits <- y_orig                  # 0s and 1s
dipper_inits <- function() list(phi = runif(1,0,1),
                                  p = runif(1,0,1),
                                  z = zinits)
head(zinits)


## -----------------------------------------------------------------------------
samples <- nimbleMCMC(dipper_code_basic, dipper_constants, dipper_data,
                      dipper_inits(), samplesAsCodaMCMC = TRUE)


## -----------------------------------------------------------------------------
plot(samples)


## -----------------------------------------------------------------------------
zinits <- y_orig
y <- y_orig
zdata <- matrix(NA, nrow = nrow(y), ncol = ncol(y))
for(i in 1:nrow(zinits)) {
  known_alive <- range(which(zinits[i,] == 1))
  zinits[i, known_alive[1] : known_alive[2] ] <- NA
  zdata[i, known_alive[1] : known_alive[2] ] <- 1
}
head(zinits)
head(zdata)


## -----------------------------------------------------------------------------
dipper_data$z <- zdata
dipper_inits <- function() list(phi = runif(1,0,1),
                                p = runif(1,0,1),
                                z = zinits)
samples <- nimbleMCMC(dipper_code_basic, dipper_constants, dipper_data,
                      dipper_inits(), samplesAsCodaMCMC = TRUE)
# This will run a bit faster.


## ----dipper_code_dcat---------------------------------------------------------
dipper_code_dcat <- nimbleCode({
  phi ~ dunif(0, 1) # prior survival
  p ~ dunif(0, 1) # prior detection
  # likelihood
  gamma[1,1:2] <- c(phi, 1-phi)      # Pr(alive t -> alive t+1), Pr(alive t -> dead t+1)
  gamma[2,1:2] <- c(0, 1)            # Pr(dead t -> alive t+1), Pr(dead t -> dead t+1)
  delta[1:2] <- c(1, 0)              # Pr(alive t = 1) = 1, Pr(dead t = 1) = 0
  omega[1,1:2] <- c(1 - p, p)        # Pr(alive t -> non-detected t), Pr(alive t -> detected t)
  omega[2,1:2] <- c(1, 0)            # Pr(dead t -> non-detected t), Pr(dead t -> detected t)
  for (i in 1:N){
    z[i,first[i]] ~ dcat(delta[1:2]) # Illustrates initial state probabilities
    for (j in (first[i]+1):T){
      z[i,j] ~ dcat(gamma[z[i,j-1], 1:2])
      y[i,j] ~ dcat(omega[z[i,j], 1:2])
    }
  }
})


## -----------------------------------------------------------------------------
dipper_data$y <- y_orig
dipper_data$y <- dipper_data$y + 1 # 1 = unobserved. 2 = observed.
zinits_dcat <- zinits
zinits_dcat[ zinits_dcat == 0] <- 2 # 1 = alive. 2 = dead.
head(dipper_data$y)
head(dipper_data$z)  # inits and data for z complement each other
head(zinits_dcat)


## -----------------------------------------------------------------------------
dipper_inits <- function() list(phi = runif(1,0,1),
                                p = runif(1,0,1),
                                z = zinits_dcat)
samples <- nimbleMCMC(dipper_code_dcat, dipper_constants, dipper_data,
                      dipper_inits(), samplesAsCodaMCMC = TRUE)


## -----------------------------------------------------------------------------
plot(samples)


## ----dipper_code_dCJS---------------------------------------------------------
dipper_code_dCJS <- nimbleCode({
  phi ~ dunif(0, 1) # survival prior
  p ~ dunif(0, 1)   # detection prior
  # likelihood
  for (i in 1:N){
    y[i, first[i]:T] ~ dCJS_ss(probSurvive=phi, probCapture=p, len=T-first[i]+1)
  }
})


## -----------------------------------------------------------------------------
dipper_inits <- function() list(phi = runif(1,0,1),
                                p = runif(1,0,1))
dipper_data$z <- NULL   # No latent states.
dipper_data$y <- y_orig # Back to 0 = unobserved, 1 = observed
samples <- nimbleMCMC(dipper_code_dCJS, dipper_constants, dipper_data,
                      dipper_inits(), samplesAsCodaMCMC = TRUE)


## -----------------------------------------------------------------------------
plot(samples)


## -----------------------------------------------------------------------------
y <- c(1, 0, 1, 1, 0, 0)
probSurvive <- c(0.7, 0.6, 0.5, 0.4, 0.3)     # survival from t to t+1
probCapture <- c(NA, 0.5, 0.7, 0.5, 0.6, 0.4) # capture probability at t
dCJS_vv(y, probSurvive = probSurvive, probCapture = probCapture)
# (In a model, we need explicit indexing with [].  Outside of a model, we don't.)
c_dCJS_vv <- compileNimble(dCJS_vv)
c_dCJS_vv(y, probSurvive = probSurvive, probCapture = probCapture)


## -----------------------------------------------------------------------------
dCJS_vv


## ---- eval=FALSE--------------------------------------------------------------
## # Run this code chunk only manually, since it has debugging in it.
## Rmodel <- nimbleModel(dipper_code_dCJS, dipper_constants, dipper_data, dipper_inits())
## Rmcmc <- buildMCMC(Rmodel)
## debugonce(dCJS_ss)
## Rmcmc$run(niter = 4)

