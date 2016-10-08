
# This script takes the output from Step1_models_v3 type runs (now v8 and v9)
# then calculates RF ranks and plots them, as well as an idea of LASSO ranks.

# FM_RF_L_VI & other RF are lists of variable importance of number of repeats 
# for cv by number of folds per repeat. Access to 3rd fold of 20th repeat is FM_RF_L_VI[[20]][[3]]

# FM_LS_Q_a0 is a simple matrix of intercept values for lasso, 
# dimensions are number of repeats x number of folds 

# FM_LS_Q_beta has a similar structure to FM_RF_L_VI and the like. list of lists, nrepeats x nfolds

rm(list = ls())
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,e1071,
               psych, foreach, randomForest, doParallel, inTrees,tableone,caret,
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
load("LS_v9.rda")
load("RF_v9.rda")

# parameters
nreps <- 26
nFeat <- 54

#####################################################
####      Calculate Ranks from RF results        ####
#####################################################
# size of VI for v9 is 54x26
VIs <- as.data.frame(matrix(0,nFeat,nreps))
rownames(VIs) <- rownames(FM_RF_L_VI[[1]][[1]])
for (i in 1:nreps)
{
  gg <- sort.int(FM_RF_L_VI[[1]][[i]], decreasing = TRUE, index.return = TRUE)
  for (j in 1:nFeat)
  {
    VIs[j,i] <- which(gg$ix == j)
  }
}

## Now get the ranks
ranks <- matrix(0,nFeat,1)
for (i in 1:nFeat)
{
  ranks[i] <- exp(mean(log(as.matrix(VIs[i,])), trim = 0.05))
}
ranks <- as.data.frame(ranks, row.names = rownames(VIs))

#####################################################
####   Calculate and compare predictions LS, RF  ####
#####################################################
# It's meaningless to do this in LOO-CV, so we do in v8
load("LS_v8.rda")
load("RF_v8.rda")

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
