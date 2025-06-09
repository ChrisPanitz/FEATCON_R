createStimArrays <- function(inputData, csCoding = "feature") {
  nrCS = max(inputData$csThisTime)
  nrTrials = max(inputData$trial)
  
  startingWeights <- inputData$contRating[1:nrCS] # only makes sense for combo coding
  startingWeights <- startingWeights[order(inputData$csThisTime[1:nrCS])]
  
  inputData <- inputData[inputData$trial>0,]
  
  usArray <- inputData$us#[(nrCS+1):(nrTrials+nrCS)]
  csArray <- array(0, dim = c(nrCS, nrTrials))
  
  if (csCoding == "combo") {
    for (tr in 1:nrTrials) {
      csArray[inputData$cs1[tr],tr] <- 1
    }
    csNames <- c("C+P+", "C+P-", "C-P+", "C-P-") # currently hard-coded for 4 CS
  }
  
  if (csCoding == "feature") {
    for (tr in (nrCS+1):(nrTrials+nrCS)) {
      csArray[inputData$colFeat[tr],(tr-nrCS)] <- 1
      csArray[inputData$patFeat[tr]+2,(tr-nrCS)] <- 1
    }
    csNames <- c("Col+", "Col-", "Pat+", "Pat-") # currently hard-coded for 2 features
  }
  
  outList <- list(csArray = csArray,
                  usArray = usArray,
                  startingWeights = startingWeights,
                  csNames = csNames)
  
  return(outList)
}
