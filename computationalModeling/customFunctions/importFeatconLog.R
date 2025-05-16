importFeatconLog <- function(filename) {
  
  inputData <- read.table(file = filename,
                          header = TRUE,
                          sep = ";",
                          dec = ".",
                          strip.white = TRUE)
  
  # CHECK AGAIN
  inputData$colFeat[inputData$cs1 == 1 | inputData$cs1 == 2] <- 1
  inputData$colFeat[inputData$cs1 == 3 | inputData$cs1 == 4] <- 2
  inputData$patFeat[inputData$cs1 == 1 | inputData$cs1 == 3] <- 1
  inputData$patFeat[inputData$cs1 == 2 | inputData$cs1 == 4] <- 2
  
  return(inputData)

}

