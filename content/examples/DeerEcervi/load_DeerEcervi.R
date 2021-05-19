library(nimble)
DeerEcervi <- read.table('DeerEcervi.txt', header = TRUE)

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

DEconstants <- list(num_farms = 24,
                    num_animals = 826,
                    length = DeerEcervi$ctrLength,
                    sex = DeerEcervi$Sex,
                    farm_ids = DeerEcervi$farm_ids)

DEdata <- list(Ecervi_01 = DeerEcervi$Ecervi_01)

DEinits <- function() {
  list(sex_int = c(0, 0),
       length_coef = c(0, 0),
       farm_sd = 1,
       farm_effect = rnorm(24, 0, 1) )
}

set.seed(123)
DEinits_vals <- DEinits()

## JAGS-compatible version
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

