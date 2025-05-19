# load functions
library(here)
parentFolder <- here()

library(ggplot2)

source(paste0(parentFolder, "/computationalModeling/customFunctions/importFeatconLog.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/createStimArrays.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/kalmanFilter.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/optimizeParametersKalman.R"))


loadfile = paste0(parentFolder,"/computationalModeling/logfiles/105_FEATCON_logfile.txt")
partData <- importFeatconLog(loadfile)
partDataNoHab <- partData[partData$trial>0,]

stimArrays <- createStimArrays(partData, csCoding = "feature")

optiParam <- optimizeParametersKalman(empData = partDataNoHab$contRating,
                                      csArray = stimArrays$csArray,
                                      usArray = stimArrays$usArray,
                                      startW = 0,
                                      expectationTime = "post",
                                      tauSqFixed = FALSE, sigmaSqRFixed = FALSE,
                                      tauSq = 0.01, sigmaSqR = 0.25) # start values for optimization

optiParam <- optimizeParametersKalman(empData = partDataNoHab$contRating[partDataNoHab$cs1 == 4],
                                      csArray = stimArrays$csArray[,partDataNoHab$cs1 == 4],
                                      usArray = stimArrays$usArray[partDataNoHab$cs1 == 4],
                                      startW = 0,
                                      expectationTime = "post",
                                      tauSqFixed = FALSE, sigmaSqRFixed = FALSE,
                                      tauSq = 0.01, sigmaSqR = 0.25) # start values for optimization
  
optiParam$par

optiKalman <- kalmanFilter(csArray = stimArrays$csArray,
                           usArray = stimArrays$usArray,
                           startW = 0,
                           #tauSq = 0.001, sigmaSqR = 0.20)
                           tauSq = optiParam$par["tauSq"],
                           sigmaSqR = optiParam$par["sigmaSqR"])

###
# ggplot() +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[1,]), color = "red") +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[2,]), color = "green") +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[3,])) +
#   geom_point(aes(x = 0:160, y = optiKalman$assocStrength[3,]), shape = "circle") +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[4,])) +
#   geom_point(aes(x = 0:160, y = optiKalman$assocStrength[4,]), shape = "diamond")

###


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
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==1]), color = "black") +
  geom_line(aes(x = 1:40, y = optiKalman$outcomeExpPost[partDataNoHab$cs1==1]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_pp,2)))), color = "black") +
  geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiKalman$tauSq,2)))), color = "black") +
  geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiKalman$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Kalman Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSpp

plotCSpm <- ggplot() +
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==2]), color = "black") +
  geom_line(aes(x = 1:40, y = optiKalman$outcomeExpPost[partDataNoHab$cs1==2]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==2] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_pm,2)))), color = "black") +
  geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiKalman$tauSq,2)))), color = "black") +
  geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiKalman$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Kalman Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSpm

plotCSmp <- ggplot() +
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==3]), color = "black") +
  geom_line(aes(x = 1:40, y = optiKalman$outcomeExpPost[partDataNoHab$cs1==3]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==3] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_mp,2)))), color = "black") +
  geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiKalman$tauSq,2)))), color = "black") +
  geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiKalman$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Kalman Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSmp

plotCSmm <- ggplot() +
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==4]), color = "black") +
  geom_line(aes(x = 1:40, y = optiKalman$outcomeExpPost[partDataNoHab$cs1==4]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==4] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_mm,2)))), color = "black") +
  geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiKalman$tauSq,2)))), color = "black") +
  geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiKalman$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Kalman Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSmm

