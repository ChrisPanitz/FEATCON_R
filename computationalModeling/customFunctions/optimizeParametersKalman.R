optimizeParametersKalman <- function(empData, csArray, usArray, expectationTime, tauSqFixed = FALSE, sigmaSqRFixed = FALSE, tauSq, sigmaSqR, startW = NULL) {

  if (tauSqFixed + sigmaSqRFixed == 2) {
    stop("Only fixed parameters - no optimization possible. Set tauSqFixed, sigmaSqRFixed, or both to FALSE.")
  }
  
  fitKalmanWrapper <- function(params, csArray, usArray, empData) {
    if (tauSqFixed + sigmaSqRFixed == 0) {
      # Unpack parameters
      tauSq <- params[1]
      sigmaSqR <- params[2]
  
      # Run Kalman filter
      model <- kalmanFilter(csArray = csArray, 
                            usArray = usArray,
                            startW = startW,
                            tauSq = tauSq, 
                            sigmaSqR = sigmaSqR)
    
    } else if (sigmaSqRFixed == TRUE) {
      # Unpack parameters
      tauSq <- params[1]
  
      # Run Kalman filter
      model <- kalmanFilter(csArray = csArray, 
                            usArray = usArray, 
                            startW = startW,
                            tauSq = tauSq,
                            sigmaSqR = sigmaSqR)
    
    } else if (tauSqFixed == TRUE) {
      # Unpack parameters
      sigmaSqR <- params[1]
        
      # Run Kalman filter
      model <- kalmanFilter(csArray = csArray, 
                            usArray = usArray, 
                            startW = startW,
                            tauSq = tauSq,
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
  if (tauSqFixed + sigmaSqRFixed == 0) {
    startParams <- c(tauSq = tauSq, sigmaSqR = sigmaSqR)
    
  } else if (sigmaSqRFixed == TRUE) {
    startParams <- c(tauSq = tauSq)
    
    
  } else if (tauSqFixed == TRUE) {
    startParams <- c(sigmaSqR = sigmaSqR)
    
  }

  # Fit model
  fit <- optim(par = startParams, 
               fn = fitKalmanWrapper, 
               csArray = csArray, 
               usArray = usArray, 
               empData = empData, 
               method = "L-BFGS-B",
               lower = c(1e-6, 1e-6),  # prevent negative diffusion parameter and variances
               upper = c(1, 0.5))
  
  return(fit)
}
# Extract best-fitting parameters
#fit$par

#checkResults <- kalmanFilter(csArray = csArray, usArray = usArray, tauSq = fit$par["tauSq"], sigmaSqR = fit$par["sigmaSqR"])

#cor.test(x = checkResults$outcomeExpPost, y = empData)
#ggplot() +
#  geom_line(aes(x = 1:length(checkResults$outcomeExpPost), y = scale(checkResults$outcomeExpPost)), color = "red") +
#  geom_line(aes(x = 1:length(empData), y = scale(empData)), color = "green") 
