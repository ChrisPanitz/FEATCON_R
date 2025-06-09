### KALMAN FILTER SCRIPT ###
# Written based on formulas in Atlas et al (2022, Psychophysiology) and Li et al. (2011, Nature Neuroscience)
hybridPearceHall <- function(csArray, usArray, eta = 0.5, kappaScaling = 1, startAlpha = 1, startW = 0, csNames = NULL) {

  ## initializing variables
  nrCS <- dim(csArray)[1]
  nrTrials <- dim(csArray)[2]

  # w = associative strength
  w <- array(data = NA, dim = c(nrCS,nrTrials+1))
  
  # initial associative strengths set to 0
  if (is.null(startW)) {
    startW <- rep(0,nrCS)
  }
  if (is.null(csNames)) {
    csNames <- c()
    for (i in 1:nrCS) {csNames[i] <- paste0("cs",i)}
  }
  w[,1] <- startW

  # empty arrays for...
  # alpha (dynamic associability)
  alphaArray <- array(NA, dim = c(nrCS,nrTrials+1))
  alphaArray[,1] <- startAlpha
  
  # outcome expectation (combined for all present CS)
  v <- array(NA, dim = nrTrials)
  vPost <- array(NA, dim = nrTrials)
  
  # prediction error
  delta <- array(NA, dim = nrTrials)

  for (i in 1:nrTrials) {
    # compute outcome expectation 
    v[i] <- matrix(w[,i],1) %*% matrix(csArray[,i])

    # compute delta
    delta[i] <- usArray[i] - v[i]
  
    # update alpha
    alphaArray[,i+1] <- alphaArray[,i]
    alphaArray[which(csArray[,i]==1),i+1] <- eta*abs(delta[i]) + (1-eta)*alphaArray[which(csArray[,i]==1),i]

    # update associative strength 
    w[,i+1] <- w[,i]
    w[which(csArray[,i]==1),i+1] <- w[which(csArray[,i]==1),i] + kappaScaling*alphaArray[which(csArray[,i]==1),i]*delta[i]
    
    # compute updated outcome expectation
    vPost[i] <- matrix(w[,i+1],1) %*% matrix(csArray[,i])
  }
  
  outList <- list(csNames = csNames,
                  trial = 1:nrTrials+1,
                  csArray = csArray,
                  usArray = usArray,
                  assocStrength = w,
                  alpha = alphaArray,
                  outcomeExp = v,
                  outcomeExpPost = vPost,
                  predError = delta,
                  kappa = kappaScaling,
                  eta = eta,
                  startAlpha = startAlpha,
                  startAssocStrenght = startW)
  
  return(outList)
}