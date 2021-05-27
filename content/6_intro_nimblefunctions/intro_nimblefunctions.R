## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(nimble)


## -----------------------------------------------------------------------------
calcDeps <- nimbleFunction(
  setup = function(model, nodes) { # setup function gives first stage of evalution
    calcNodes <- model$getDependencies(nodes)
  },
  run = function() {               # run function (or other methods) give second stage of evaluation
    ans <- model$calculate(calcNodes)
    return(ans)
    returnType(double())
  }
)


## -----------------------------------------------------------------------------
regrCode <- nimbleCode({
  b0 ~ dnorm(0, sd = 100)
  b1 ~ dnorm(0, sd = 100)
  sigma ~ dunif(0, 100)
  for(i in 1:n)
     y[i] ~ dnorm(b0 + b1*x[i], sd = sigma)
})

set.seed(1)
n <- 100
x <- runif(n)
b0_true <- 0.3
b1_true <- 0.5
sigma_true <- 0.25
y <- rnorm(n, b0_true + b1_true*x, sigma_true)

regrModel <- nimbleModel(regrCode, data = list(y = y), constants = list(n = n, x = x),
                     inits = list(b0 = 1, b1 = 0, sigma = 1))


## -----------------------------------------------------------------------------
calcDeps_regr <- calcDeps(regrModel, c('b0', 'b1', 'sigma'))
calcDeps_regr$run()   ## second stage: run code


## -----------------------------------------------------------------------------
class(calcDeps_regr) ## We could have used nimbleFunction's name argument to set the class name
calcDeps_regr$calcNodes[1:20]


## ---- eval=FALSE--------------------------------------------------------------
## calcDeps <- nimbleFunction(
##   setup = function(model, nodes) {
##     browser()
##     calcNodes <- model$getDependencies(nodes)
##   },
##   run = function() {
##     browser()
##     ans <- model$calculate(calcNodes)
##     return(ans)
##     returnType(double())
##   }
## ) ## warning about not being able to compiled with browser() is expected.


## ---- eval=FALSE--------------------------------------------------------------
## calcDeps_regr <- calcDeps(model, c('b0', 'b1', 'sigma')) ## We'll see the setup code followed by internal code.
## calcDeps_regr$run()


## ---- regr-objective, eval----------------------------------------------------
objective <- nimbleFunction(
    setup = function(model, nodes) {
        calcNodes <- model$getDependencies(nodes)
        elements <- model$expandNodeNames(nodes, returnScalarComponents = TRUE)
        n <- length(elements)
    },
    run = function(par = double(1)) {
        returnType(double(0))
        if(length(par) != n)
           stop("Input length does not match number of parameter elements.")
        values(model, nodes) <<- par   # assignment into non-local (nf) variables 
        ans <- model$calculate(calcNodes)  # local assignment
        return(ans)
    }
)


## ---- regr-specialized--------------------------------------------------------
rObjective <- objective(regrModel, c('b0', 'b1', 'sigma'))
cRegrModel <- compileNimble(regrModel)   # remember to compile model first
cObjective <- compileNimble(rObjective, project = regrModel)


## ---- regr-optimize-----------------------------------------------------------
set.seed(1)
system.time(optR <- optim(c(0, 0, 1), rObjective$run, control = list(fnscale = -1)))
system.time(optC <- optim(c(0, 0, 1), cObjective$run, control = list(fnscale = -1)))
optR
optC

