# This script takes the output from Step1_models_v3 type runs (v8-v13)
# then calculates RF ranks and plots them, as well as an idea of LASSO ranks.

# FM_RF_L_VI & other RF are lists of variable importance of number of repeats 
# for cv by number of folds per repeat. Access to 3rd fold of 20th repeat is FM_RF_L_VI[[20]][[3]]

# FM_LS_Q_a0 is a simple matrix of intercept values for lasso, 
# dimensions are number of repeats x number of folds 

# FM_LS_Q_beta has a similar structure to FM_RF_L_VI and the like. list of lists, nrepeats x nfolds

# rm(list = ls()) #clear the workspace
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
load("RF_v8.rda")

# parameters
nReps <- length(FM_RF_L_VI)
nFolds <- length(FM_RF_L_VI[[1]])
nFeat <- length(FM_RF_L_VI[[1]][[1]])

#####################################################
####      Calculate Ranks from RF results        ####
#####################################################
# first let's reorder all the variable importance vectors in a feature x repetition matrix
# remember: repetition = nReps x nFolds for repeated nFold cross validation

# Initialize VI matrix
VIs_FM <- as.data.frame(matrix(0,nFeat,nReps*nFolds))
rownames(VIs_FM) <- rownames(FM_RF_L_VI[[1]][[1]])

VIs_WO <- as.data.frame(matrix(0,nFeat,nReps*nFolds))
rownames(VIs_WO) <- rownames(WO_RF_L_VI[[1]][[1]])

for (i in 1:nReps)
{
  for (kk in 1:nFolds)
  {
    coln <- (nFolds*(i-1))+kk # column number
    VIs_FM[,coln] <- FM_RF_L_VI[[i]][[kk]]
  }
}

for (i in 1:nReps)
{
  for (kk in 1:nFolds)
  {
    coln <- (nFolds*(i-1))+kk # column number
    VIs_WO[,coln] <- WO_RF_L_VI[[i]][[kk]]
  }
}

write.csv(VIs_FM,file = 'RF_VIs_FM_Raw_v8.csv')
write.csv(VIs_WO,file = 'RF_VIs_WO_Raw_v8.csv')
