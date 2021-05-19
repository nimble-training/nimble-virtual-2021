## custom MCMC sampler
sampler_conditional_RW <- nimbleFunction(
    name = 'sampler_conditional_RW',
    contains = sampler_BASE,    
    setup = function(model, mvSaved, target, control) {
        ## control list extraction
        logScale            <- extractControlElement(control, 'log',                 FALSE)
        adaptive            <- extractControlElement(control, 'adaptive',            TRUE)
        adaptInterval       <- extractControlElement(control, 'adaptInterval',       200)
        adaptFactorExponent <- extractControlElement(control, 'adaptFactorExponent', 0.8)
        scale               <- extractControlElement(control, 'scale',               1)
        ## set incorrect defaults so that it throws an error if not explicitly set
        index               <- extractControlElement(control, 'index',               1:2)
        
        ## checks on target
        if(length(target) != 1) stop('length of target not equal to one')
        targetAsScalar <- model$expandNodeNames(target, returnScalarComponents = TRUE)
        d <- length(targetAsScalar)
        if(d <= 1) stop('target does not have more than one dimension')
        
        ## checks on index
        if(length(index) != 1) stop('length of index must be 1')
        if(index < 1 | index > d) stop('index must be within number of dimensions of target')
        
        ## node list generation
        calcNodes <- model$getDependencies(target)
        calcNodesNoSelf <- model$getDependencies(target, self = FALSE)
        isStochCalcNodesNoSelf <- model$isStoch(calcNodesNoSelf)   ## should be made faster
        calcNodesNoSelfDeterm <- calcNodesNoSelf[!isStochCalcNodesNoSelf]
        calcNodesNoSelfStoch <- calcNodesNoSelf[isStochCalcNodesNoSelf]
        ## numeric value generation
        scaleOriginal <- scale
        timesRan      <- 0
        timesAccepted <- 0
        timesAdapted  <- 0
        scaleHistory  <- c(0, 0)   ## scaleHistory
        acceptanceHistory  <- c(0, 0)   ## scaleHistory
        if(nimbleOptions('MCMCsaveHistory')) {
            saveMCMChistory <- TRUE
        } else saveMCMChistory <- FALSE
        optimalAR     <- 0.44
        gamma1        <- 0
        ## checks
        if(model$isDiscrete(target))     stop('cannot use RW sampler on discrete-valued target; try slice sampler')
        if(adaptFactorExponent < 0)      stop('cannot use RW sampler with adaptFactorExponent control parameter less than 0')
        if(scale < 0)                    stop('cannot use RW sampler with scale control parameter less than 0')
    },
    run = function() {
        currentValue <- model[[target]]
        propValue <- model[[target]]
        propLogScale <- 0
        if(logScale) { 
            propLogScale <- rnorm(1, mean = 0, sd = scale)
            propValue[index] <- currentValue[index] * exp(propLogScale)
        } else {
            propValue[index] <- rnorm(1, mean = currentValue[index],  sd = scale)
        }
        model[[target]] <<- propValue
        logMHR <- model$calculateDiff(target)
        if(logMHR == -Inf) {
            nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
            ## Drawing a random number is needed during first testing
            ## of this step in order to keep the random numbers identical
            ## to old behavior to see if tests that depend on particular
            ## sample sequences pass.  Rather than calling runif(1, 0, 1) here,
            ## we call decide() to ensure same behavior.
            ## jump <- decide(logMHR)
            ## When new behavior is acceptable, we can remove the above line
            ## and uncomment the following:
            jump <- FALSE
        } else {
            logMHR <- logMHR + model$calculateDiff(calcNodesNoSelf) + propLogScale
            jump <- decide(logMHR)
            if(jump) {
                nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
                nimCopy(from = model, to = mvSaved, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
            } else {
                nimCopy(from = mvSaved, to = model, row = 1, nodes = target, logProb = TRUE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfDeterm, logProb = FALSE)
                nimCopy(from = mvSaved, to = model, row = 1, nodes = calcNodesNoSelfStoch, logProbOnly = TRUE)
            }
        }
        if(adaptive)     adaptiveProcedure(jump)
    },
    methods = list(
        adaptiveProcedure = function(jump = logical()) {
            timesRan <<- timesRan + 1
            if(jump)     timesAccepted <<- timesAccepted + 1
            if(timesRan %% adaptInterval == 0) {
                acceptanceRate <- timesAccepted / timesRan
                timesAdapted <<- timesAdapted + 1
                if(saveMCMChistory) {
                    setSize(scaleHistory, timesAdapted)                 ## scaleHistory
                    scaleHistory[timesAdapted] <<- scale                ## scaleHistory
                    setSize(acceptanceHistory, timesAdapted)            ## scaleHistory
                    acceptanceHistory[timesAdapted] <<- acceptanceRate  ## scaleHistory
                }
                gamma1 <<- 1/((timesAdapted + 3)^adaptFactorExponent)
                gamma2 <- 10 * gamma1
                adaptFactor <- exp(gamma2 * (acceptanceRate - optimalAR))
                scale <<- scale * adaptFactor
                timesRan <<- 0
                timesAccepted <<- 0
            }
        },
        getScaleHistory = function() {       ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(scaleHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
                return(numeric(1, 0))
            }
        },          
        getAcceptanceHistory = function() {  ## scaleHistory
            returnType(double(1))
            if(saveMCMChistory) {
                return(acceptanceHistory)
            } else {
                print("Please set 'nimbleOptions(MCMCsaveHistory = TRUE)' before building the MCMC")
                return(numeric(1, 0))
            }
        },
        reset = function() {
            scale <<- scaleOriginal
            timesRan      <<- 0
            timesAccepted <<- 0
            timesAdapted  <<- 0
            if(saveMCMChistory) {
                scaleHistory  <<- c(0, 0)    ## scaleHistory
                acceptanceHistory  <<- c(0, 0)
            }
            gamma1 <<- 0
        }
    )
)

