optimizeParametersPearceHall <- function(empData, csArray, usArray, expectationTime,
                                         optimCriterion = "corr", maxIterations = 200,
                                         etaFixed = FALSE, kappaFixed = FALSE, startAlphaFixed = FALSE, startWFixed = FALSE,
                                         eta = 0.5, kappaScaling = 1, startAlpha = 1, startW = 0,
                                         etaBound = c(1e-6,1), kappaBound = c(1e-6,1), startAlphaBound = c(-1,1), startWBound = c(-1,1)) {

  # only run if there is at least one free parameter
  if (etaFixed + kappaFixed + startAlphaFixed + startWFixed == 4) {
    stop("Only fixed parameters - no optimization possible. Set at least one of etaFixed, kappaFixed, stratAlphaFixed, or startWFixed to FALSE.")
  }
  
  # load DEoptim package
  library(DEoptim)
  
  # check what parameters are free and list parameter names + boundaries
  paramNames <- NULL
  lowerBounds <- NULL
  upperBounds <- NULL
  
  if (etaFixed == FALSE) {
    paramNames <- c(paramNames, "eta")
    lowerBounds <- c(lowerBounds, etaBound[1])
    upperBounds <- c(upperBounds, etaBound[2])}
  if (kappaFixed == FALSE) {
    paramNames <- c(paramNames, "kappa")
    lowerBounds <- c(lowerBounds, kappaBound[1])
    upperBounds <- c(upperBounds, kappaBound[2])}
  if (startAlphaFixed == FALSE) {
    paramNames <- c(paramNames, "startAlpha")
    lowerBounds <- c(lowerBounds, startAlphaBound[1])
    upperBounds <- c(upperBounds, startAlphaBound[2])}
  if (startWFixed == FALSE) {
    paramNames <- c(paramNames, "startW")
    lowerBounds <- c(lowerBounds, startWBound[1])
    upperBounds <- c(upperBounds, startWBound[2])}
  
  # define wrapper function to run DEoptim on hybrid Pearce-Hall model function
  fitPearceHallWrapper <- function(params, csArray, usArray, empData) {
    nrPar <- 1
    if (etaFixed == FALSE) {eta <- params[nrPar]; nrPar <- nrPar+1}
    if (kappaFixed == FALSE) {kappa <- params[nrPar]; nrPar <- nrPar+1}
    if (startAlphaFixed == FALSE) {startAlpha <- params[nrPar]; nrPar <- nrPar+1}
    if (startWFixed == FALSE) {startW <- params[nrPar]}
      

    # Run hybrid Pearce-Hall model
    model <- hybridPearceHall(csArray = csArray,
                              usArray = usArray,
                              startW = startW,
                              eta = eta,
                              kappaScaling = kappaScaling,
                              startAlpha = startAlpha)
    
    # Model predictions (outcome expectations)
    if (expectationTime == "pre") {
      predictions <- model$outcomeExp
    }
    else if (expectationTime == "post") {
      predictions <- model$outcomeExpPost
    }
    
    # Return what metric to optimize?
    if (optimCriterion == "rmse") {
      # Error metric: RMSE
      rmse <- sqrt(mean((empData - predictions)^2))
      return(rmse)      
    }
    else if (optimCriterion == "corr") { 
      # "Error" metric: r
      corr <- cor(empData,predictions)
      return(1-abs(corr))
    }
  }

  # Fit model
  fit <- DEoptim(fn = fitPearceHallWrapper,
                 csArray = csArray,
                 usArray = usArray,
                 empData = empData,
                 lower = lowerBounds,
                 upper = upperBounds,
                 control = list(itermax = maxIterations))

  
  names(fit$optim$bestmem) <- paramNames
  return(fit)
}