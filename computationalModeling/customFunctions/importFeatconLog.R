importFeatconLog <- function(filename) {
  
  inputData <- read.table(file = filename,
                          header = TRUE,
                          sep = ";",
                          dec = ".",
                          strip.white = TRUE)
  
  # CHECK AGAIN
  inputData$colFeat[inputData$csBl1 == 1 | inputData$csBl1 == 2] <- 1
  inputData$colFeat[inputData$csBl1 == 3 | inputData$csBl1 == 4] <- 2
  inputData$patFeat[inputData$csBl1 == 1 | inputData$csBl1 == 3] <- 1
  inputData$patFeat[inputData$csBl1 == 2 | inputData$csBl1 == 4] <- 2
  
  return(inputData)

}

