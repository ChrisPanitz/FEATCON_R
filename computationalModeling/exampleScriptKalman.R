# load functions
library(here)
parentFolder <- here()

library(ggplot2)

source(paste0(parentFolder, "/computationalModeling/customFunctions/importFeatconLog.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/createStimArrays.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/kalmanFilter.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/optimizeParametersKalman.R"))


loadfile = paste0(parentFolder,"/computationalModeling/logfiles/103_FEATCON_logfile.txt")
partData <- importFeatconLog(loadfile)
partDataNoHab <- partData[partData$trial>0,]

stimArrays <- createStimArrays(partData, csCoding = "combo")

optiParam <- optimizeParametersKalman(empData = partDataNoHab$contRating,
                                      csArray = stimArrays$csArray,
                                      usArray = stimArrays$usArray,
                                      startW = 0,
                                      expectationTime = "post",
                                      tauSqFixed = FALSE, sigmaSqRFixed = FALSE,
                                      tauSq = 0.01, sigmaSqR = 0.20) # start values for optimization
  
optiParam$par

optiKalman <- kalmanFilter(csArray = stimArrays$csArray,
                           usArray = stimArrays$usArray,
                           startW = 0,
                           tauSq = optiParam$par["tauSq"],
                           sigmaSqR = optiParam$par["sigmaSqR"])

rmse <- sqrt(mean((partDataNoHab$contRating - optiKalman$outcomeExpPost)^2))
corr <- cor(partDataNoHab$contRating, optiKalman$outcomeExpPost)
corr_pp <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 1], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 1])
corr_pm <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 2], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 2])
corr_mp <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 3], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 3])
corr_mm <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 4], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 4])



# plot all trials
# plotAllCS <- ggplot() +
#                geom_line(aes(x = 1:160, y = scale(partDataNoHab$contRating)), color = "black") +
#                geom_line(aes(x = 1:160, y = scale(optiKalman$outcomeExpPost)), color = "blue")
# plotAllCS


# plot all trials
plotCSpp <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==1])), color = "black") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==1])), color = "blue") +
  #geom_vline(xintercept = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), color = "gray") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), y = 1.5)) +
  geom_label(aes(x = 40, y = 1.3, label = paste0("r = ", as.character(round(corr_pp,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "z-values of predicted (blue) and observed (black) responses") +
  theme_classic()
plotCSpp

plotCSpm <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==2])), color = "black") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==2])), color = "blue") +
  #geom_vline(xintercept = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), color = "gray") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), y = 1.5)) +
  geom_label(aes(x = 40, y = 1.3, label = paste0("r = ", as.character(round(corr_pm,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "z-values of predicted (blue) and observed (black) responses") +
  theme_classic()
plotCSpm

plotCSmp <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==3])), color = "black") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==3])), color = "blue") +
  #geom_vline(xintercept = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), color = "gray") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==3] == 1), y = 1.5)) +
  geom_label(aes(x = 40, y = 1.3, label = paste0("r = ", as.character(round(corr_mp,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "z-values of predicted (blue) and observed (black) responses") +
  theme_classic()
plotCSmp

plotCSmm <- ggplot() +
  geom_line(aes(x = 1:40, y = scale(partDataNoHab$contRating[partDataNoHab$cs1==4])), color = "black") +
  geom_line(aes(x = 1:40, y = scale(optiKalman$outcomeExpPost[partDataNoHab$cs1==4])), color = "blue") +
  #geom_vline(xintercept = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), color = "gray") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==4] == 1), y = 1.5)) +
  geom_label(aes(x = 40, y = 1.3, label = paste0("r = ", as.character(round(corr_mm,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "z-values of predicted (blue) and observed (black) responses") +
  theme_classic()
plotCSmm
