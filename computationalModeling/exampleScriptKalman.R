# load functions
library(here)
parentFolder <- here()

library(ggplot2)

source(paste0(parentFolder, "/computationalModeling/customFunctions/importFeatconLog.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/createStimArrays.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/kalmanFilter.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/optimizeParametersKalman.R"))

partVec = c(101,103:111)
partDF = data.frame(partID = integer(length = length(partVec)),
                    tauSq = double(length = length(partVec)),
                    sigmaSqR = double(length = length(partVec)),
                    rmse = double(length = length(partVec)),
                    corr = double(length = length(partVec)))

for (partI in 1:length(partVec)) {
  loadfile = paste0(parentFolder,"/computationalModeling/logfiles/",as.character(partVec[partI]), "_FEATCON_logfile.txt")
  partData <- importFeatconLog(loadfile)
  partDataNoHab <- partData[partData$trial>0,]
  
  stimArrays <- createStimArrays(partData, csCoding = "combo")
  
  optiParam <- optimizeParametersKalman(empData = partDataNoHab$contRating,
                                        csArray = stimArrays$csArray,
                                        usArray = stimArrays$usArray,
                                        startW = 0,
                                        expectationTime = "post",
                                        tauSqFixed = FALSE, sigmaSqRFixed = TRUE,
                                        tauSq = 0.01, sigmaSqR = 0.20) # start values for optimization
  
  # optiParam <- optimizeParametersKalman(empData = partDataNoHab$contRating[partDataNoHab$cs1 == 4],
  #                                       csArray = stimArrays$csArray[,partDataNoHab$cs1 == 4],
  #                                       usArray = stimArrays$usArray[partDataNoHab$cs1 == 4],
  #                                       startW = 0,
  #                                       expectationTime = "post",
  #                                       tauSqFixed = FALSE, sigmaSqRFixed = TRUE,
  #                                       tauSq = 0.01, sigmaSqR = 0.25) # start values for optimization
    
  optiParam$par
  
  optiKalman <- kalmanFilter(csArray = stimArrays$csArray,
                             usArray = stimArrays$usArray,
                             startW = 0,
                             #tauSq = 0.001, sigmaSqR = 0.20)
                             tauSq = optiParam$par["tauSq"],
                             sigmaSqR = 0.20) #optiParam$par["sigmaSqR"])
  
  ### model fit parameters ###
  nrObs <- length(optiKalman$outcomeExpPost)
  residuals <- partDataNoHab$contRating - optiKalman$outcomeExpPost
  nrFreeParam <- length(optiParam$par)
  # root of the mean squared error
  rmse <- sqrt(mean((residuals)^2))
  # Pearson correlation
  corr <- cor(partDataNoHab$contRating, optiKalman$outcomeExpPost)
  # log-likelihood
  residVar <- var(residuals) # residual variance
  log_likelihood <- -0.5 * nrObs * log(2 * pi * residVar) - sum(residuals^2) / (2 * residVar)
  # Bayesian Information Criterion (BIC)
  BIC <- nrFreeParam * log(nrObs) - 2*log_likelihood
  
  
  partDF$partID[partI] <- partVec[partI]
  partDF$tauSq[partI] <- optiParam$par["tauSq"]
  partDF$sigmaSqR[partI] <- optiParam$par["sigmaSqR"]
  partDF$rmse[partI] <- rmse
  partDF$corr[partI] <- round(corr,3)
  
  
} # end participant loop


###
# ggplot() +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[1,]), color = "red") +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[2,]), color = "green") +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[3,])) +
#   geom_point(aes(x = 0:160, y = optiKalman$assocStrength[3,]), shape = "circle") +
#   geom_line(aes(x = 0:160, y = optiKalman$assocStrength[4,])) +
#   geom_point(aes(x = 0:160, y = optiKalman$assocStrength[4,]), shape = "diamond")

###



corr_pp <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 1], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 1])
corr_pm <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 2], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 2])
corr_mp <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 3], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 3])
corr_mm <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 4], optiKalman$outcomeExpPost[partDataNoHab$cs1 == 4])



# plot all trials
plotAllCS <- ggplot() +
               geom_line(aes(x = 1:160, y = scale(partDataNoHab$contRating)), color = "black") +
               geom_line(aes(x = 1:160, y = scale(optiKalman$outcomeExpPost)), color = "blue")
plotAllCS


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

