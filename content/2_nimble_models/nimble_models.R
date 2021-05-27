## ----setup, include=FALSE-----------------------------------------------------
library(nimble)
source(file.path("..", "examples", "DeerEcervi", "load_DeerEcervi.R"), 
       chdir = TRUE)


## -----------------------------------------------------------------------------
DEcode_partial <- nimbleCode({
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
  # <snip>
  })


## ---- eval=FALSE--------------------------------------------------------------
## farm_sd ~ dunif(0, 20)
## farm_prec <- 1/farm_sd^2
## for(i in 1:num_farms) {
##   farm_effect[i] ~ dnorm(0, farm_prec)
## }


## ---- eval=FALSE--------------------------------------------------------------
## for(i in 1:num_animals) {
##   logit(disease_probability[i]) <-
##     sex_int[ sex[i] ] +
##     length_coef[ sex[i] ] * cLength[i] +
##     farm_effect[ farm_ids[i] ]
##   Ecervi_01[i] ~ dbern(disease_probability[i])
## }


## ---- eval=FALSE--------------------------------------------------------------
## for(j in 1:num_farms)
##   for(i in 1:num_animals_per_farm) {
##     logit(disease_probability[j, i]) <- farm_effect[j]
##     Ecervi_01[j, i] ~ dbern(disease_probability[j, i])
##   }


## ----eval = FALSE-------------------------------------------------------------
## x[1:5] <- (A[1:5, 1:5] %*% b[1:5] + c[1:5])


## ----eval = FALSE-------------------------------------------------------------
## x[1:5] <- (A[1:5, 1:5] %*% b[1:5] + c[1:5])[,1]


## ----eval=FALSE---------------------------------------------------------------
## code <- nimbleCode({
##   sigma ~ dunif(0, 10)
##   beta0 ~ dnorm(0, sd = 1000)
##   beta1 ~ dnorm(0, sd = 1000)
##   if(INCLUDE_X2) {
##     beta2 ~ dnorm(0, sd = 1000)
##   }
##   for(i in 1:10) {
##     if(INCLUDE_X2) {
##       y[i] ~ dnorm(beta0 + beta1 * x1[i] + beta2 * x2[i], sd = sigma)
##     } else {
##       y[i] ~ dnorm(beta0 + beta1 * x1[i], sd = sigma)
##     }
##   }
## })
## 
## INCLUDE_X2 <- FALSE
## m1 <- nimbleModel(code)
## INCLUDE_X2 <- TRUE
## m2 <- nimbleModel(code)


## -----------------------------------------------------------------------------
set.seed(1)
code <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 100)
  for(i in 1:4) {
    predicted.y[i] <- intercept + slope * x[i]
    y[i] ~ dnorm(predicted.y[i], sd = sigma)
  }
})
model <- nimbleModel(code, 
                     data = list(y = rnorm(4)),
                     inits = list(intercept = 0.5, 
                                  slope = 0.2, 
                                  sigma = 1,
                                  x = c(0.1, 0.2, 0.3, 0.4)))


## ---- linmodel-graph, echo = FALSE--------------------------------------------
layout <- matrix(ncol = 2, byrow = TRUE,
   # These seem to be rescaled to fit in the plot area,
   # so I'll just use 0-100 as the scale
                 data = c(33, 100,
                          66, 100,
                          50, 0, # first three are parameters
                          15, 50, 35, 50, 55, 50, 75, 50, # x's
                          20, 75, 40, 75, 60, 75, 80, 75, # predicted.y's
                          25, 25, 45, 25, 65, 25, 85, 25) # y's
                 )

sizes <- c(45, 30, 30,
           rep(20, 4),
           rep(50, 4),
           rep(20, 4))

edge.color <- "black"
stoch.color <- "deepskyblue2"
det.color <- "orchid3"
rhs.color <- "gray73"
fill.color <- c(
    rep(stoch.color, 3),
    rep(rhs.color, 4),
    rep(det.color, 4),
    rep(stoch.color, 4)
)

plot(model$graph, vertex.shape = "crectangle",
     vertex.size = sizes,
     vertex.size2 = 20,
     layout = layout,
     vertex.label.cex = 1.0,
     vertex.color = fill.color,
     edge.width = 3,
     asp = 0.5,
     edge.color = edge.color)


## -----------------------------------------------------------------------------
model$sigma
model$x
model$x[3] <- 0.6
model$x


## -----------------------------------------------------------------------------
model$y
model$y[1] <- 0.8
model$y


## -----------------------------------------------------------------------------
model$getNodeNames()


## -----------------------------------------------------------------------------
model$getNodeNames(dataOnly = TRUE)


## -----------------------------------------------------------------------------
model$getNodeNames(determOnly = TRUE)


## -----------------------------------------------------------------------------
model$isData('y')
model$isData('x')


## -----------------------------------------------------------------------------
model$getDependencies("x[2]")


## -----------------------------------------------------------------------------
model$getDependencies("sigma")


## -----------------------------------------------------------------------------
model$getDependencies("slope")


## ---- eval=FALSE--------------------------------------------------------------
## for(i in 1:num_animals) {
##   logit(disease_probability[i]) <-
##     sex_int[ sex[i] ] +
##     length_coef[ sex[i] ] * Length[i] +
##     farm_effect[ farm_ids[i] ]
##   Ecervi_01[i] ~ dbern(disease_probability[i])


## -----------------------------------------------------------------------------
DEconstants <- list(num_farms = 24,
                    num_animals = 826,
                    cLength = DeerEcervi$cLength,
                    sex = DeerEcervi$Sex,
                    farm_ids = DeerEcervi$farm_ids)
DEdata <- list(Ecervi_01 = DeerEcervi$Ecervi_01)

m_light <- nimbleModel(DEcode, constants = DEconstants)

DEconstants <- list(num_farms = 24,
                    num_animals = 826,
                    cLength = DeerEcervi$cLength)
DEdata <- list(Ecervi_01 = DeerEcervi$Ecervi_01,
                    sex = DeerEcervi$Sex,
                    farm_ids = DeerEcervi$farm_ids)
m_heavy <-  nimbleModel(DEcode, data = DEdata, constants = DEconstants)


## -----------------------------------------------------------------------------
m_light$getDependencies('farm_effect[1]')
too_many_deps <- m_heavy$getDependencies('farm_effect[1]')
length(too_many_deps)
head(too_many_deps, n = 30)


## -----------------------------------------------------------------------------
code2 <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma ~ dunif(0, 100)
  predicted.y[1:4] <- intercept + slope * x[1:4] # vectorized node
  for(i in 1:4) {
    # predicted.y[i] <- intercept + slope * x[i] # scalar nodes (earlier model version)
    y[i] ~ dnorm(predicted.y[i], sd = sigma)
  }
})

model2 <- nimbleModel(code2, 
                      data = list(y = rnorm(4)),
                      inits = list(intercept = 0.5, 
                                   slope = 0.2, 
                                   sigma = 1,
                                   x = c(0.1, 0.2, 0.3, 0.4)))


## -----------------------------------------------------------------------------
model2$getNodeNames()


## -----------------------------------------------------------------------------
model2$getDependencies('x[2]')


## -----------------------------------------------------------------------------
model2$calculate('y[1:4]')


## -----------------------------------------------------------------------------
model2$getDependencies('intercept')
model2$calculate(model2$getDependencies('intercept'))


## -----------------------------------------------------------------------------
model2$sigma
model2$simulate('sigma')
model2$sigma


## -----------------------------------------------------------------------------
model2$y
model2$simulate('y') ## Will not over-write data nodes
model2$y
model2$simulate('y', includeData = TRUE) ## will over-write data nodes
model2$y


## -----------------------------------------------------------------------------
code3 <- nimbleCode({
  intercept ~ dnorm(0, sd = 1000)
  slope ~ dnorm(0, sd = 1000)
  sigma2 ~ dinvgamma(1, 1)  # this sort of prior not generally recommended
  for(i in 1:4) {
    y[i] ~ dnorm(intercept + slope * x[i], var = sigma2)
  }
})
model3 <- nimbleModel(code3, 
                      data = list(y = rnorm(4)),
                      inits = list(intercept = 0.5, 
                                   slope = 0.2, 
                                   sigma2 = 1,
                                   x = c(0.1, 0.2, 0.3, 0.4)))


## -----------------------------------------------------------------------------
model3$getNodeNames()


## ---- lifted------------------------------------------------------------------
model3$sigma2 <- 100
model3$lifted_sqrt_oPsigma2_cP
model3$simulate('y', includeData = TRUE)
summary(model3$y)
depNodes <- model3$getDependencies('sigma2', self = FALSE)
depNodes
model3$simulate(depNodes, includeData = TRUE)
model3$lifted_sqrt_oPsigma2_cP
summary(model3$y)


## ---- generic-simulate--------------------------------------------------------
simulate_downstream <- function(model, nodes) {
  downstream_nodes <- model$getDependencies(nodes, downstream = TRUE)
  model$simulate( downstream_nodes, includeData = TRUE )
  logProb <- model$calculate( downstream_nodes )
  logProb
}


## -----------------------------------------------------------------------------
model3$y
simulate_downstream(model3, 'sigma2')
model3$y


## ---- eval=FALSE--------------------------------------------------------------
## model$farm_sd <- 1
## model$sex_int <- c(1.5, 2.5)
## model$length_coef <- c(0.3, 0.5)
## # explore the arguments for `getDependencies` and `simulate`
## deps <- model$getDependencies('???')
## model$simulate(???)

