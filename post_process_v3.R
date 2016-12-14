
# This script takes the output from Step1_models_v3 type runs (v8-v13)
# then calculates RF ranks and plots them, as well as an idea of LASSO ranks.

# FM_RF_L_VI & other RF are lists of variable importance of number of repeats 
# for cv by number of folds per repeat. Access to 3rd fold of 20th repeat is FM_RF_L_VI[[20]][[3]]

# FM_LS_Q_a0 is a simple matrix of intercept values for lasso, 
# dimensions are number of repeats x number of folds 

# FM_LS_Q_beta has a similar structure to FM_RF_L_VI and the like. list of lists, nrepeats x nfolds

rm(list = ls()) #clear the workspace
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,e1071,combinat,
               psych, foreach, randomForest, doParallel, inTrees,tableone,caret,mgcv,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)
# options(error=recover)

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations

#####################################################
####      Load the results from Step1_models     ####
#####################################################

source("functions.R")
source("create_df1_df2_no_unknown.R")

# Since Saria recommended LOO-CV, we'll take v9 for now
# load("2016dec12_RF_nReps100_nFolds13.rda")

# parameters
# nReps <- length(FM_RF_L_VI)
# nFolds <- length(FM_RF_L_VI[[1]])
# nFeat <- length(FM_RF_L_VI[[1]][[1]])
# 
# #####################################################
# ####      Calculate Ranks from RF results        ####
# #####################################################
# # first let's reorder all the variable importance vectors in a feature x repetition matrix
# # remember: repetition = nReps x nFolds for repeated nFold cross validation
# 
# # Initialize VI matrix
# VIs_FM <- as.data.frame(matrix(0,nFeat,nReps*nFolds))
# rownames(VIs_FM) <- rownames(FM_RF_L_VI[[1]][[1]])
# 
# VIs_WO <- as.data.frame(matrix(0,nFeat,nReps*nFolds))
# rownames(VIs_WO) <- rownames(WO_RF_L_VI[[1]][[1]])
# 
# for (i in 1:nReps)
# {
#   for (kk in 1:nFolds)
#   {
#     coln <- (nFolds*(i-1))+kk # column number
#     gg <- sort.int(FM_RF_L_VI[[i]][[kk]], decreasing = TRUE, index.return = TRUE)
#     for (j in 1:nFeat)
#     {
#       VIs_FM[j,coln] <- which(gg$ix == j)
#     }
#   }
# }
# 
# for (i in 1:nReps)
# {
#   for (kk in 1:nFolds)
#   {
#     coln <- (nFolds*(i-1))+kk # column number
#     gg <- sort.int(WO_RF_L_VI[[i]][[kk]], decreasing = TRUE, index.return = TRUE)
#     for (j in 1:nFeat)
#     {
#       VIs_WO[j,coln] <- which(gg$ix == j)
#     }
#   }
# }
# 
# ## Now get the ranks
# ranks_FM <- matrix(0,nFeat,1)
# ranks_WO <- matrix(0,nFeat,1)
# 
# for (i in 1:nFeat)
# {
#   ranks_FM[i] <- exp(mean(log(as.matrix(VIs_FM[i,])), trim = 0.05))
#   ranks_WO[i] <- exp(mean(log(as.matrix(VIs_WO[i,])), trim = 0.05))
# }
# ranks_FM <- as.data.frame(ranks_FM, row.names = rownames(VIs_FM))
# ranks_WO <- as.data.frame(ranks_WO, row.names = rownames(VIs_WO))
# 
# save(ranks_FM,ranks_WO,VIs_FM,VIs_WO, file = "RF_Ranks_nReps100_nFolds13.rda")
#####################################################
####   Calculate and compare predictions LS, RF  ####
#####################################################
# It's meaningless to do this in LOO-CV, so we do in v8
load("2016dec12_LS_nReps1_nFolds26.rda")
load("2016dec12_RF_nReps1_nFolds26.rda")

nreps <- dim(FM_LS_L)[2]
RMSE_FM_LS_L <- matrix(0,nreps,1)
RMSE_FM_LS_Q <- matrix(0,nreps,1)
RMSE_FM_RF_L <- matrix(0,nreps,1)
RMSE_FM_RF_Q <- matrix(0,nreps,1)

RMSE_WO_LS_L <- matrix(0,nreps,1)
RMSE_WO_LS_Q <- matrix(0,nreps,1)
RMSE_WO_RF_L <- matrix(0,nreps,1)
RMSE_WO_RF_Q <- matrix(0,nreps,1)

# 100 predictions each, matrix format

for (i in 1:nreps)
{
  FM_LS_L_mod <- lm(yFM~FM_LS_L[,i])
  RMSE_FM_LS_L[i] <- round(sqrt(mean(resid(FM_LS_L_mod)^2)), 2)

  FM_LS_Q_mod <- lm(yFM~FM_LS_Q[,i])
  RMSE_FM_LS_Q[i] <- round(sqrt(mean(resid(FM_LS_Q_mod)^2)), 2)

  FM_RF_L_mod <- lm(yFM~FM_RF_L[,i])
  RMSE_FM_RF_L[i] <- round(sqrt(mean(resid(FM_RF_L_mod)^2)), 2)

  FM_RF_Q_mod <- lm(yFM~FM_RF_Q[,i])
  RMSE_FM_RF_Q[i] <- round(sqrt(mean(resid(FM_RF_Q_mod)^2)), 2)

  # Wolf
  WO_LS_L_mod <- lm(yWO~WO_LS_L[,i])
  RMSE_WO_LS_L[i] <- round(sqrt(mean(resid(WO_LS_L_mod)^2)), 2)

  WO_LS_Q_mod <- lm(yWO~WO_LS_Q[,i])
  RMSE_WO_LS_Q[i] <- round(sqrt(mean(resid(WO_LS_Q_mod)^2)), 2)

  WO_RF_L_mod <- lm(yWO~WO_RF_L[,i])
  RMSE_WO_RF_L[i] <- round(sqrt(mean(resid(WO_RF_L_mod)^2)), 2)

  WO_RF_Q_mod <- lm(yWO~WO_RF_Q[,i])
  RMSE_WO_RF_Q[i] <- round(sqrt(mean(resid(WO_RF_Q_mod)^2)), 2)
}

heights <- matrix(0,8,1)
heights[1] <- mean(RMSE_FM_LS_L)
heights[2] <- mean(RMSE_FM_LS_Q)
heights[3] <- mean(RMSE_FM_RF_L)
heights[4] <- mean(RMSE_FM_RF_Q)

heights[5] <- mean(RMSE_FM_LS_L)
heights[6] <- mean(RMSE_FM_LS_Q)
heights[7] <- mean(RMSE_FM_RF_L)
heights[8] <- mean(RMSE_FM_RF_Q)

heights <- as.data.frame(heights)

rownames(heights) <- c("Fugl-Meyer lasso linear", "Fugl-Meyer lasso quad","Fugl-Meyer RF linear",
                       "Fugl-Meyer RF quad","Wolf lasso linear","Wolf lasso quad","Wolf RF linear",
                       "Wolf RF quad")
barplot(as.matrix(t(heights)))
