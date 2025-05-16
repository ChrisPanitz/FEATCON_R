# load functions
source("customFunctions/importFeatconLog.R")
source("customFunctions/createStimArrays.R")
source("customFunctions/kalmanFilter.R")
source("customFunctions/optimizeParametersKalman.R")
library(here)
library(ggplot2)


parentFolder <- here()
loadfile = paste0(parentFolder,"/computationalModeling/logfiles/101_FEATCON_logfile.txt")
partData <- importFeatconLog(loadfile)
partDataNoHab <- partData[partData$trial>0,]

stimArrays <- createStimArrays(partData, csCoding = "feature")

optiParam <- optimizeParametersKalman(empData = partDataNoHab$contRating,
                                      csArray = stimArrays$csArray,
                                      usArray = stimArrays$usArray,
                                      expectationTime = "post",
                                      startTau = 0.01, startSigma = 0.20)
  
optiParam$par

optiKalman <- kalmanFilter(csArray = stimArrays$csArray,
                           usArray = stimArrays$usArray,
                           tauSq = optiParam$par["tauSq"],
                           sigmaSqR = optiParam$par["sigmaSqR"])

rmse <- sqrt(mean((partDataNoHab$contRating - optiKalman$outcomeExpPost)^2))
corr <- cor(partDataNoHab$contRating, optiKalman$outcomeExpPost)



# plot all trials
plotAllCS <- ggplot() +
               geom_line(aes(x = 1:160, y = scale(partDataNoHab$contRating)), color = "darkred") +
               geom_line(aes(x = 1:160, y = scale(optiKalman$outcomeExpPost)), color = "darkgreen")
plotAllCS


# plot all trials
plotCSpp <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==1])), color = "darkred") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==1])), color = "darkgreen")
plotCSpp

plotCSpm <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==2])), color = "darkred") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==2])), color = "darkgreen")
plotCSpm

plotCSmp <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==3])), color = "darkred") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==3])), color = "darkgreen")
plotCSmp

plotCSmm <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==4])), color = "darkred") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==4])), color = "darkgreen")
plotCSmm
