## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,
                      cache = TRUE)
has_ggplot2 <- require(ggplot2)
has_mcmcplots <- require(mcmcplots)
has_coda <- require(coda)
generate_original_results <- FALSE


## -----------------------------------------------------------------------------
library(nimble)


## -----------------------------------------------------------------------------
DeerEcervi <- read.table(file.path('..', 'examples', 'DeerEcervi', 'DeerEcervi.txt'), header = TRUE)
summary(DeerEcervi)

## Create presence/absence data from counts.
DeerEcervi$Ecervi_01 <- DeerEcervi$Ecervi
DeerEcervi$Ecervi_01[DeerEcervi$Ecervi>0] <- 1
## Set up naming convention for centered and uncentered lengths for exercises later
DeerEcervi$unctrLength <- DeerEcervi$Length
## Center Length for better interpretation
DeerEcervi$ctrLength <- DeerEcervi$Length - mean(DeerEcervi$Length)
## Make a factor version of Sex for plotting
DeerEcervi$fSex <- factor(DeerEcervi$Sex)
## Make a factor and id version of Farm
DeerEcervi$fFarm <- factor(DeerEcervi$Farm)
DeerEcervi$farm_ids <- as.numeric(DeerEcervi$fFarm)




## -----------------------------------------------------------------------------
DEcode <- nimbleCode({
  for(i in 1:2) {
    # Priors for intercepts and length coefficients for sex = 1 (male), 2 (female)
    sex_int[i] ~ dnorm(0, sd = 1000)
    length_coef[i] ~ dnorm(0, sd = 1000)
  }

  # Priors for farm random effects and their standard deviation.
  farm_sd ~ dunif(0, 20)
  for(i in 1:num_farms) {
    farm_effect[i] ~ dnorm(0, sd = farm_sd)
  }

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
DEconstants <- list(num_farms = 24,
                    num_animals = 826,
                    length = DeerEcervi$ctrLength,
                    sex = DeerEcervi$Sex,
                    farm_ids = DeerEcervi$farm_ids)

DEmodel <- nimbleModel(DEcode,
                       constants = DEconstants)


## -----------------------------------------------------------------------------
DEmodel$setData(list(Ecervi_01 = DeerEcervi$Ecervi_01))
# This sets the values and *flags the nodes as data*.
DEinits <- function() {
  list(sex_int = c(0, 0),
       length_coef = c(0, 0),
       farm_sd = 1,
       farm_effect = rnorm(24, 0, 1) )
}

set.seed(123)
DEmodel$setInits(DEinits())


## -----------------------------------------------------------------------------
DEmcmc <- buildMCMC(DEmodel)


## -----------------------------------------------------------------------------
cDEmodel <- compileNimble(DEmodel) 
# First call to compileNimble in a session is slower than later calls.
cDEmcmc <- compileNimble(DEmcmc, project = DEmodel)


## -----------------------------------------------------------------------------
cDEmcmc$run(10000)


## -----------------------------------------------------------------------------
samples1 <- as.matrix(cDEmcmc$mvSamples)


## ----eval=FALSE---------------------------------------------------------------
## # Run this code if you want to generate your own results.
## # They won't over-write results that come with these slides.
## library(mcmcplots)
## mcmcplot(samples1, dir = ".", filename = "Ecervi_samples_mcmcplot")




## ----eval = FALSE-------------------------------------------------------------
## # We haven't provided coda figures, but you can make make them if you want.
## library(coda)
## pdf("Ecervi_samples_coda.pdf")
## plot(as.mcmc(samples1))
## dev.off()


## -----------------------------------------------------------------------------
set.seed(123)
DEdataAndConstants <- c(DEconstants, 
                        list(Ecervi_01 = DeerEcervi$Ecervi_01))
samples2 <- nimbleMCMC(DEcode,
                       constants = DEdataAndConstants,
                       inits = DEinits,
                       niter = 10000,
                       nburnin = 1000,
                       nchains = 2,
                       samplesAsCodaMCMC = TRUE)
summary(samples2) ## from coda


## -----------------------------------------------------------------------------
samples3 <- runMCMC(cDEmcmc, 
                    niter = 10000,
                    nburnin = 1000,
                    nchains = 2,
                    samplesAsCodaMCMC = TRUE)
summary(samples3)

