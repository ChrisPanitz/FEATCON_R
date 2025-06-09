optimizeParametersKalman <- function(empData, csArray, usArray, expectationTime,
                                     optimCriterion = "corr", maxIterations = 200,
                                     tauSqFixed = FALSE, sigmaSqRFixed = FALSE, startWFixed = FALSE,
                                     tauSq = 0.01, sigmaSqR = 1, startW = 0,
                                     tauSqBound = c(1e-6,1), sigmsSqRBound = c(1e-6,1), startWBound = c(1e-6,1)) {

  # only run if there is at least one free parameter
  if (tauSqFixed + sigmaSqRFixed + startWFixed == 3) {
    stop("Only fixed parameters - no optimization possible. Set at least on of tauSqFixed, sigmaSqRFixed, and startWFixed to FALSE.")
  }
  
  # load DEoptim package
  library(DEoptim)
  
  # check what parameters are free and list parameter names + boundaries
  paramNames <- NULL
  lowerBounds <- NULL
  upperBounds <- NULL
  
  if (tauSqFixed == FALSE) {
    paramNames <- c(paramNames, "tauSq")
    lowerBounds <- c(lowerBounds, tauSqBound[1])
    upperBounds <- c(upperBounds, tauSqBound[2])}
  if (sigmaSqRFixed == FALSE) {
    paramNames <- c(paramNames, "sigmaSqR")
    lowerBounds <- c(lowerBounds, sigmaSqRBound[1])
    upperBounds <- c(upperBounds, sigmaSqRBound[2])}
  if (startWFixed == FALSE) {
    paramNames <- c(paramNames, "startW")
    lowerBounds <- c(lowerBounds, startWBound[1])
    upperBounds <- c(upperBounds, startWBound[2])}
  
  # define wrapper function to run DEoptim on Kalman filter function
  fitKalmanWrapper <- function(params, csArray, usArray, empData) {
    nrPar <- 1
    if (tauSqFixed == FALSE) {tauSq <- params[nrPar]; nrPar <- nrPar+1}
    if (sigmaSqRFixed == FALSE) {sigmaSqR <- params[nrPar]; nrPar <- nrPar+1}
    if (startWFixed == FALSE) {startW <- params[nrPar]}
      

    # Run Kalman filter
    model <- kalmanFilter(csArray = csArray,
                          usArray = usArray,
                          startW = startW,
                          tauSq = tauSq,
                          sigmaSqR = sigmaSqR)
    
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
  fit <- DEoptim(fn = fitKalmanWrapper,
                 csArray = csArray,
                 usArray = usArray,
                 empData = empData,
                 lower = lowerBounds,
                 upper = upperBounds,
                 control = list(itermax = maxIterations))

  
  names(fit$optim$bestmem) <- paramNames
  return(fit)
}