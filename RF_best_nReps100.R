## This function 

#########################
## Setup and libraries ##
#########################

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
## Data Folder: Dropbox/Research/Import            ##
##                 Data Prep - Importing           ##
#####################################################

source("functions.R")
load('2016dec12_Basics.rda')
nReps <- 100
nFolds <- 4
nSubj <- 26
nFeat <- 51       # Number of features of linear models
nFeat2 <- 1365    # Number of features of quadratic models

# save the predictions
FM_RF_L <- matrix(0,nSubj,nReps)
FM_RF_Q <- matrix(0,nSubj,nReps)
WO_RF_L <- matrix(0,nSubj,nReps)
WO_RF_Q <- matrix(0,nSubj,nReps)

# save the feature importance ranks
FM_RF_L_VI <- matrix(0,nFeat ,nReps)
FM_RF_Q_VI <- matrix(0,nFeat2,nReps)
WO_RF_L_VI <- matrix(0,nFeat ,nReps)
WO_RF_Q_VI <- matrix(0,nFeat2,nReps)

# Numbers in this loop
# 1: FM L || 2: FM Q || 3: RF L || 4: RF Q

for (i in 1:nReps)
{
  print(i)
  # Build the forests
  RF_model1 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% randomForest(df1,yFM,ntree=ntree, corr.bias=FALSE)
  RF_model2 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% randomForest(df3,yFM,ntree=ntree, corr.bias=FALSE)
  RF_model3 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% randomForest(df1,yWO,ntree=ntree, corr.bias=FALSE)
  RF_model4 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% randomForest(df3,yWO,ntree=ntree, corr.bias=FALSE)
  print(paste("Created forests for repeat ",i))
  
  # Calculate variable importance
  FM_RF_L_VI[,i] <- RF_model1$importance
  FM_RF_Q_VI[,i] <- RF_model2$importance
  WO_RF_L_VI[,i] <- RF_model3$importance
  WO_RF_Q_VI[,i] <- RF_model4$importance
  print(paste("Calculated Variable Importance for repeat ",i))
  
  # Calculate predictions for each run
  FM_RF_L[,i] <- predict(RF_model1, df1, type = "response" , predict.all=FALSE)
  FM_RF_Q[,i] <- predict(RF_model2, df3, type = "response" , predict.all=FALSE)
  WO_RF_L[,i] <- predict(RF_model3, df1, type = "response" , predict.all=FALSE)
  WO_RF_Q[,i] <- predict(RF_model4, df3, type = "response" , predict.all=FALSE)
  print(paste("Calculated Predictions for repeat ",i))
}

write.csv(FM_RF_L,"2016dec12_FM_RF_L_Best_nReps100.csv")
write.csv(FM_RF_Q,"2016dec12_FM_RF_Q_Best_nReps100.csv")
write.csv(WO_RF_L,"2016dec12_WO_RF_L_Best_nReps100.csv")
write.csv(WO_RF_Q,"2016dec12_WO_RF_Q_Best_nReps100.csv")

save(FM_RF_L,FM_RF_Q,WO_RF_L,WO_RF_Q,FM_RF_L_VI,FM_RF_Q_VI,WO_RF_L_VI,WO_RF_Q_VI, file = "2016dec12_RF_Best_nReps100.rda")