optimizeParametersKalman <- function(empdata, expectationTime) {

  fitKalmanWrapper <- function(params, csArray, usArray, empData) {
    # Unpack parameters
    tauSq <- params[1]
    sigmaSqR <- params[2]

    # Run Kalman filter
    model <- kalmanFilter(csArray = csArray, 
                          usArray = usArray, 
                          tauSq = tauSq, 
                          sigmaSqR = sigmaSqR)
    
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
  startParams <- c(tauSq = 0.01, sigmaSqR = 0.2)
  
  # Fit model
  fit <- optim(par = startParams, 
               fn = fitKalmanWrapper, 
               csArray = csArray, 
               usArray = usArray, 
               empData = empData, 
               method = "L-BFGS-B",
               lower = c(1e-6, 1e-6),  # prevent negative diffusion parameter and variances
               upper = c(1, 10))
}
# Extract best-fitting parameters
#fit$par

#checkResults <- kalmanFilter(csArray = csArray, usArray = usArray, tauSq = fit$par["tauSq"], sigmaSqR = fit$par["sigmaSqR"])

#cor.test(x = checkResults$outcomeExpPost, y = empData)
#ggplot() +
#  geom_line(aes(x = 1:length(checkResults$outcomeExpPost), y = scale(checkResults$outcomeExpPost)), color = "red") +
#  geom_line(aes(x = 1:length(empData), y = scale(empData)), color = "green") 
