# This version adds bootstrapping to LASSO, so the model predictions
# can be averaged and mean and std can be calculated for these
# predictions, the results will be forwarded to matlab for plotting

#########################
## Setup and libraries ##
#########################

rm(list = ls())
# options(error=recover)
library(MASS)
library(R.matlab)
library(corrplot)
library(GGally)
library(tableone)
library(lars)
library(glmnet)
library(coefplot)
library(ggplot2)
library(foreach)
library(randomForest)
library(doParallel)
library(inTrees)
library(qpcR)

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations

#####################################################
## Data Folder: Dropbox/Research/Import            ##
##                 Data Prep - Importing           ##
#####################################################

source("functions.R")
source("create_df1_df2.R")

######################################################
## From this point forward, the numbers attached to ##
## variables will mean the following:               ##
##  # 1: Linear model, Fugl-Meyer                   ##
##  # 2: Linear Model, partial Fugl-Meyer           ##
##  # 3: Linear Model, Wolf                         ##
##  # 4: Quadratic Model, Fugl-Meyer                ##
##  # 5: Quadratic Model, partial Fugl-Meyer        ##
##  # 6: Quadratic Model, Wolf                      ##
######################################################

#################################################
## 05 - LASSO Cross-validation & Bootstrapping ##
#################################################
## Mod 1: linear
## Mod 2: quadratic
predLASSO1<-NULL
predLASSO2<-NULL
predLASSO3<-NULL
predLASSO4<-NULL
predLASSO5<-NULL
predLASSO6<-NULL

for (i in 1:100)
{
  predLASSO1[[i]] <- cv_mod(feat=Features,out=yFM,stand=1,isWO = 0)
  predLASSO2[[i]] <- cv_mod(feat=Features,out=yPartFM,stand=1,isWO = 0)
  predLASSO3[[i]] <- cv_mod(feat=Features,out=yWO,stand=1,isWO = 1)
  predLASSO4[[i]] <- cv_mod(feat=df2,out=yFM,stand=1,isWO = 0)
  predLASSO5[[i]] <- cv_mod(feat=df2,out=yPartFM,stand=1,isWO = 0)
  predLASSO6[[i]] <- cv_mod(feat=df2,out=yWO,stand=1,isWO = 1)
}


# LASSO_pred1 <- cleanup_LASSO(predLASSO1,1)

LASSO_pred1 <- cleanup_LASSO(predLASSO1)
LASSO_pred2 <- cleanup_LASSO(predLASSO2)
LASSO_pred3 <- cleanup_LASSO(predLASSO3)
LASSO_pred4 <- cleanup_LASSO(predLASSO4)
LASSO_pred5 <- cleanup_LASSO(predLASSO5)
LASSO_pred6 <- cleanup_LASSO(predLASSO6)

save(LASSO_pred1,LASSO_pred2,LASSO_pred3,LASSO_pred4,LASSO_pred5,LASSO_pred6,file = "lasso_pred_corrbias.rda")

#######################################
## 06 - Random Forests Bootstrapping ##
#######################################
predRF1 <- NULL
predRF2 <- NULL
predRF3 <- NULL
predRF4 <- NULL
predRF5 <- NULL
predRF6 <- NULL

for (i in 1:100)
{
  predRF1[[i]] <- cv_RF(feat=Features, out=yFM)
  predRF2[[i]] <- cv_RF(feat=Features, out=yPartFM)
  predRF3[[i]] <- cv_RF(feat=Features, out=yWO)
  predRF4[[i]] <- cv_RF(feat=df2, out=yFM)
  predRF5[[i]] <- cv_RF(feat=df2, out=yPartFM)
  predRF6[[i]] <- cv_RF(feat=df2, out=yWO)
}

RF_pred1 <- cleanup_RF(predRF1)
RF_pred2 <- cleanup_RF(predRF2)
RF_pred3 <- cleanup_RF(predRF3)
RF_pred4 <- cleanup_RF(predRF4)
RF_pred5 <- cleanup_RF(predRF5)
RF_pred6 <- cleanup_RF(predRF6)

save(RF_pred1, RF_pred2, RF_pred3, RF_pred4, RF_pred5, RF_pred6, file = "RF_pred_corrbias.rda")
