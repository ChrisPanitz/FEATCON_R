### KALMAN FILTER SCRIPT ###
# Written based on formulas in Gershman (2015, PLoS Computational Biology)
kalmanFilter <- function(csArray, usArray, startW = NULL, tauSq = .01, sigmaSqR = 0.20, sigmaSqW = 1, csNames = NULL) {
  # defaults taken from Gershman (2015)
  
  ## initializing variables
  nrCS <- dim(csArray)[1]
  nrTrials <- dim(csArray)[2]

  # w = associative strength
  w <- array(data = NA, dim = c(nrCS,nrTrials+1))
  # covMat = variances and covariances of CS
  covMat <- array(data = NA, dim = c(nrCS,nrCS,nrTrials+1))

  # initial associative strengths set to 0
  if (is.null(startW)) {
    startW <- rep(0,nrCS)
  }
  if (is.null(csNames)) {
    csNames <- c()
    for (i in 1:nrCS) {csNames[i] <- paste0("cs",i)}
  }
  w[,1] <- startW

  # initial variances set to input, covariances set to 0
  covMat[, , 1] <- sigmaSqW * diag(nrCS);

  # empty arrays for...
  # Kalman gains (each CS has their own)
  kGain <- array(NA, dim = c(nrCS,nrTrials+1))
  # outcome expectation (combined for all present CS)
  v <- array(NA, dim = nrTrials)
  vPost <- array(NA, dim = nrTrials)
  # prediction error
  delta <- array(NA, dim = nrTrials+1)

  for (i in 1:nrTrials) {
    # compute Kalman gains (Gershman: Formula (11))
    kGain[,i] <- ((covMat[, , i] + tauSq*diag(nrCS)) %*% matrix(csArray[,i])) / 
    as.numeric(matrix(csArray[,i],1) %*% (covMat[, , i] + tauSq*diag(nrCS)) %*% matrix(csArray[,i]) + sigmaSqR)

    # update covariance matrix (Gershman: Formula (10))
    covMat[, , i+1] <- covMat[, , i] + tauSq*diag(nrCS) - matrix(kGain[,i]) %*% matrix(csArray[,i],1) %*% (covMat[, , i] + tauSq*diag(nrCS))
      
    # compute outcome expectation (Gershman: Formula (2))
    v[i] <- matrix(w[,i],1) %*% matrix(csArray[,i])

    # compute delta (Gershman: Page 3, text after Formula (2))
    delta[i] <- usArray[i] - v[i]

    # update associative strength (Gershman: Formula (9))
    w[,i+1] <- w[,i] + kGain[,i] * delta[i]
    
    # compute updated outcome expectation
    vPost[i] <- matrix(w[,i+1],1) %*% matrix(csArray[,i])
  }
  
  outList <- list(csNames = csNames,
                  trial = 0:nrTrials,
                  csArray = csArray,
                  usArray = usArray,
                  assocStrength = w,
                  covMat = covMat,
                  outcomeExp = v,
                  outcomeExpPost = vPost,
                  predError = delta,
                  kalmanGain = kGain,
                  tauSq = tauSq,
                  sigmaSqR = sigmaSqR)
  
  return(outList)
}