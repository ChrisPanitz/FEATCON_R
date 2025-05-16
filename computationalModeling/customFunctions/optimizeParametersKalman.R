optimizeParametersKalman <- function(empData, csArray, usArray, expectationTime, startTau, startSigma) {

  fitKalmanWrapper <- function(params, csArray, usArray, empData) {
    if (is.null(startTau)+is.null(startSigma) == 0) {
      # Unpack parameters
      tauSq <- params[1]
      sigmaSqR <- params[2]
  
      # Run Kalman filter
      model <- kalmanFilter(csArray = csArray, 
                            usArray = usArray, 
                            tauSq = tauSq, 
                            sigmaSqR = sigmaSqR)
    
    } else if (is.null(startSigma)) {
      # Unpack parameters
      tauSq <- params[1]
  
      # Run Kalman filter
      model <- kalmanFilter(csArray = csArray, 
                            usArray = usArray, 
                            tauSq = tauSq)
    
    } else if (is.null(startTau)) {
      # Unpack parameters
      sigmaSqR <- params[1]
        
      # Run Kalman filter
      model <- kalmanFilter(csArray = csArray, 
                            usArray = usArray, 
                            sigmaSqR = sigmaSqR)
      
    }
    
    # Model predictions (outcome expectations)
    if (expectationTime == "pre") {
      predictions <- model$outcomeExp
    }
    else if (expectationTime == "post") {
      predictions <- model$outcomeExpPost
    }
    
    # Error metric: RMSE
    rmse <- sqrt(mean((empData - predictions)^2))
    return(rmse)
  }

  # Starting guesses for parameters
  if (is.null(startTau)+is.null(startSigma) == 0) {
    startParams <- c(tauSq = startTau, sigmaSqR = startSigma)
    
  } else if (is.null(startSigma)) {
    startParams <- c(tauSq = startTau)
    
    
  } else if (is.null(startTau)) {
    startParams <- c(sigmaSqR = startSigma)
    
  }

  # Fit model
  fit <- optim(par = startParams, 
               fn = fitKalmanWrapper, 
               csArray = csArray, 
               usArray = usArray, 
               empData = empData, 
               method = "L-BFGS-B",
               lower = c(1e-6, 1e-6),  # prevent negative diffusion parameter and variances
               upper = c(1, 10))
  
  return(fit)
}
# Extract best-fitting parameters
#fit$par

#checkResults <- kalmanFilter(csArray = csArray, usArray = usArray, tauSq = fit$par["tauSq"], sigmaSqR = fit$par["sigmaSqR"])

#cor.test(x = checkResults$outcomeExpPost, y = empData)
#ggplot() +
#  geom_line(aes(x = 1:length(checkResults$outcomeExpPost), y = scale(checkResults$outcomeExpPost)), color = "red") +
#  geom_line(aes(x = 1:length(empData), y = scale(empData)), color = "green") 
