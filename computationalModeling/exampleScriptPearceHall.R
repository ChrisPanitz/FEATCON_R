# load functions
library(here)
parentFolder <- here()

library(ggplot2)

source(paste0(parentFolder, "/computationalModeling/customFunctions/importFeatconLog.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/createStimArrays.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/hybridPearceHall.R"))
source(paste0(parentFolder, "/computationalModeling/customFunctions/optimizeParametersPearceHall.R"))

partVec = c(101,103:111)
partDF = data.frame(partID = integer(length = length(partVec)),
                    eta = double(length = length(partVec)),
                    kappa = double(length = length(partVec)),
                    startAlpha = double(length = length(partVec)),
                    startW = double(length = length(partVec)),
                    rmse = double(length = length(partVec)),
                    corr = double(length = length(partVec)),
                    BIC = double(length = length(partVec)))

for (partI in 1:length(partVec)) {
  loadfile = paste0(parentFolder,"/computationalModeling/logfiles/",as.character(partVec[partI]), "_FEATCON_logfile.txt")
  partData <- importFeatconLog(loadfile)
  partDataNoHab <- partData[partData$trial>0,]
  
  stimArrays <- createStimArrays(partData, csCoding = "combo")
  
  optiParam <- optimizeParametersPearceHall(empData = partDataNoHab$contRating,
                                            csArray = stimArrays$csArray,
                                            usArray = stimArrays$usArray,
                                            expectationTime = "post",
                                            optimCriterion = "corr",
                                            maxIterations = 200,
                                            etaFixed = FALSE, kappaFixed = FALSE, startAlphaFixed = FALSE, startWFixed = FALSE, 
                                            etaBound = c(1e-6,1), kappaBound = c(1e-6,1), startAlphaBound = c(1e-6,1), startWBound = c(1e-6,1),
                                            eta = 0.5, kappaScaling = 1, startAlpha = 1, startW = 0)
  
  optiParam$optim$bestmem
  
  optiPearceHall <- hybridPearceHall(csArray = stimArrays$csArray,
                                      usArray = stimArrays$usArray,
                                      eta = optiParam$optim$bestmem["eta"],
                                      kappa = optiParam$optim$bestmem["kappa"],
                                      startAlpha = optiParam$optim$bestmem["startAlpha"],
                                      startW = optiParam$optim$bestmem["startW"])
  
  ### model fit parameters ###
  nrObs <- length(optiPearceHall$outcomeExpPost)
  residuals <- partDataNoHab$contRating - optiPearceHall$outcomeExpPost
  nrFreeParam <- length(optiParam$par)
  # root of the mean squared error
  rmse <- sqrt(mean((residuals)^2))
  # Pearson correlation
  corr <- cor(partDataNoHab$contRating, optiPearceHall$outcomeExpPost)
  # log-likelihood
  residVar <- var(residuals) # residual variance
  log_likelihood <- -0.5 * nrObs * log(2 * pi * residVar) - sum(residuals^2) / (2 * residVar)
  # Bayesian Information Criterion (BIC)
  BIC <- nrFreeParam * log(nrObs) - 2*log_likelihood
  
  
  partDF$partID[partI] <- partVec[partI]
  partDF$eta[partI] <- optiParam$optim$bestmem["eta"]
  partDF$kappa[partI] <- optiParam$optim$bestmem["kappa"]
  partDF$startAlpha[partI] <- optiParam$optim$bestmem["startAlpha"]
  partDF$startW[partI] <- optiParam$optim$bestmem["startW"]
  partDF$rmse[partI] <- rmse
  partDF$corr[partI] <- round(corr,3)
  partDF$BIC[partI] <- round(BIC,1)
  
  
} # end participant loop


#optiPearceHall <- hybridPearceHall(csArray = stimArrays$csArray, usArray = stimArrays$usArray, eta = -2, kappa = 22, startAlpha = -1.5, startW = 10)

###
# ggplot() +
#   geom_line(aes(x = 0:160, y = optiPearceHall$assocStrength[1,]), color = "red") +
#   geom_line(aes(x = 0:160, y = optiPearceHall$assocStrength[2,]), color = "green") +
#   geom_line(aes(x = 0:160, y = optiPearceHall$assocStrength[3,])) +
#   geom_point(aes(x = 0:160, y = optiPearceHall$assocStrength[3,]), shape = "circle") +
#   geom_line(aes(x = 0:160, y = optiPearceHall$assocStrength[4,])) +
#   geom_point(aes(x = 0:160, y = optiPearceHall$assocStrength[4,]), shape = "diamond")

###
ggplot() +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==1]), color = "red") +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==2])) +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==3])) +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==4]))

deltaP <- ggplot() + geom_line(aes(x = 1:40, y = optiPearceHall$predError[partDataNoHab$cs1==1]), color = "red")
absDeltaP <- ggplot() + geom_line(aes(x = 1:40, y = abs(optiPearceHall$predError[partDataNoHab$cs1==1])), color = "red")
alphaP <- ggplot() + geom_line(aes(x = 1:40, y = optiPearceHall$alpha[1,partDataNoHab$cs1==1]), color = "red")
wP <- ggplot() + geom_line(aes(x = 1:40, y = optiPearceHall$assocStrength[1,partDataNoHab$cs1==1]), color = "red")
vP <- ggplot() + geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==1]), color = "red")
usP <- ggplot() + geom_line(aes(x = 1:40, y = optiPearceHall$usArray[partDataNoHab$cs1==1]), color = "red")

library(ggpubr)
ggarrange(deltaP, absDeltaP, alphaP, wP, vP, usP,
          nrow = 6, ncol = 1,
          labels = c("delta","absDeltaP","alpha","w", "vP", "us"))


corr_pp <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 1], optiPearceHall$outcomeExpPost[partDataNoHab$cs1 == 1])
corr_pm <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 2], optiPearceHall$outcomeExpPost[partDataNoHab$cs1 == 2])
corr_mp <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 3], optiPearceHall$outcomeExpPost[partDataNoHab$cs1 == 3])
corr_mm <- cor(partDataNoHab$contRating[partDataNoHab$cs1 == 4], optiPearceHall$outcomeExpPost[partDataNoHab$cs1 == 4])



# plot all trials
plotAllCS <- ggplot() +
               geom_line(aes(x = 1:160, y = scale(partDataNoHab$contRating)), color = "black") +
               geom_line(aes(x = 1:160, y = scale(optiPearceHall$outcomeExpPost)), color = "blue")
plotAllCS


# plot all trials
plotCSpp <- ggplot() +
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==1]), color = "black") +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==1]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==1] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_pp,2)))), color = "black") +
  #geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiPearceHall$tauSq,2)))), color = "black") +
  #geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiPearceHall$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Pearce-Hall Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSpp

plotCSpm <- ggplot() +
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==2]), color = "black") +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==2]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==2] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_pm,2)))), color = "black") +
  #geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiPearceHall$tauSq,2)))), color = "black") +
  #geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiPearceHall$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Pearce-Hall Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSpm

plotCSmp <- ggplot() +
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==3]), color = "black") +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==3]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==3] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_mp,2)))), color = "black") +
  #geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiPearceHall$tauSq,2)))), color = "black") +
  #geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiPearceHall$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Pearce-Hall Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSmp

plotCSmm <- ggplot() +
  geom_line(aes(x = 1:40, y = partDataNoHab$contRating[partDataNoHab$cs1==4]), color = "black") +
  geom_line(aes(x = 1:40, y = optiPearceHall$outcomeExpPost[partDataNoHab$cs1==4]*100), color = "blue") +
  geom_point(aes(x = which(partDataNoHab$us[partDataNoHab$cs1==4] == 1), y = 120)) +
  geom_label(aes(x = 10, y = 5, label = paste0("r = ", as.character(round(corr_mm,2)))), color = "black") +
  #geom_label(aes(x = 20, y = 5, label = paste0("tau^2 = ", as.character(round(optiPearceHall$tauSq,2)))), color = "black") +
  #geom_label(aes(x = 30, y = 5, label = paste0("sigma^2R = ", as.character(round(optiPearceHall$sigmaSqR,2)))), color = "black") +
  scale_x_continuous(name = "trial #") +
  scale_y_continuous(name = "Probability Rating (0-100%)", breaks = seq(0,100,20),
                     sec.axis = sec_axis(name = "Pearce-Hall Prediction (blue)", trans = ~./100, breaks = seq(0.0, 1.0, 0.2))) +
  theme_classic() +
  theme(axis.line.y.right = element_line(color = "blue"),
        axis.ticks.y.right = element_line(color = "blue"),
        axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"))
plotCSmm

