# This version adds bootstrapping to LASSO, so the model predictions
# can be averaged and mean and std can be calculated for these
# predictions, the results will be forwarded to matlab for plotting

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
source("create_df1_df2_no_unknown.R")
write.csv(file="FeatureSet_v2.csv",x=df1)
#################################################
#### 05 - LASSO Cross-validation & Bootstrapping ####
#################################################

## Update Sep 2016:
# after much back and forth, the current decision is to use df1 and df3
# df2 is dropped in favor of having a model that does not extrapolate
# The numbers now are:
# 1: Linear FM LASSO
# 2: quad df3 FM LASSO
# 5: Linear WO lasso
# 6: Quad df3 WO lasso
# same for RF instead of lasso
# 
# predLASSO1      <- NULL
# predLASSO2    <-NULL
# predLASSO5      <-NULL
# predLASSO6      <-NULL
# # predLASSO5    <-NULL
# # predLASSO6      <-NULL
# # predLASSO7      <- NULL
# # predLASSO8      <- NULL
# 
# for (i in 1:100)
# {
#   predLASSO1[[i]] <- cv_mod(feat=df1,out=yFM,stand=1,isWO = 0, nf = 4)
#   predLASSO2[[i]] <- cv_mod(feat=df3,out=yFM,stand=1,isWO = 0, nf = 4)
#   predLASSO5[[i]] <- cv_mod(feat=df1,out=yWO,stand=1,isWO = 1, nf = 4)
#   predLASSO6[[i]] <- cv_mod(feat=df3,out=yWO,stand=1,isWO = 1, nf = 4)
#   # predLASSO5[[i]] <- cv_mod(feat=df2,out=yPartFM,stand=1,isWO = 0)
#   # predLASSO6[[i]] <- cv_mod(feat=df2,out=yWO,stand=1,isWO = 1, nf = 4)
#   # predLASSO7[[i]] <- cv_mod(feat=df3, out=yFM, stand=1,isWO = 0, nf=4)
#   # predLASSO8[[i]] <- cv_mod(feat=df3, out=yWO, stand=1,isWO = 1, nf=4)
# }
# 
# LASSO_pred1 <- cleanup_LASSO(predLASSO1)
# LASSO_pred2 <- cleanup_LASSO(predLASSO2)
# # LASSO_pred3 <- cleanup_LASSO(predLASSO3)
# # LASSO_pred4 <- cleanup_LASSO(predLASSO4)
# LASSO_pred5 <- cleanup_LASSO(predLASSO5)
# LASSO_pred6 <- cleanup_LASSO(predLASSO6)
# # LASSO_pred7 <- cleanup_LASSO(predLASSO7)
# # LASSO_pred8 <- cleanup_LASSO(predLASSO8)
# 
# save(LASSO_pred1,LASSO_pred2,LASSO_pred5,LASSO_pred6,file = "lasso_pred_v6.rda")

#######################################
#### 06 - Random Forests Bootstrapping ####
#######################################
# 3: linear FM RF
# 4: quad df3 FM RF
# 7: linear WO RF
# 8: quad df3 WO RF

FM_RF_L <- matrix(0,26,100)
FM_RF_Q <- matrix(0,26,100)
WO_RF_L <- matrix(0,26,100)
WO_RF_Q <- matrix(0,26,100)

for (i in 1:100)
{
  print(i)
  whoop <- createFolds(yFM, k=4, list=FALSE)
  for (j in 1:4) # 4-fold cross-validation
  {
    # first let's form the train and test sets
    train_L <- df1[whoop!=j,]; test_L <- df1[whoop==j,]; FMtrain <- yFM[whoop!=j]; WOtrain <- yWO[whoop!=j];
    train_Q <- df3[whoop!=j,]; test_Q <- df3[whoop==j,]; FMtest <- yFM[whoop==j]; WOtest <- yWO[whoop==j];
    RF_model3 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                         .packages = "randomForest") %dopar% randomForest(train_L,FMtrain,ntree=ntree, corr.bias=FALSE)
    RF_model4 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                         .packages = "randomForest") %dopar% randomForest(train_Q,FMtrain,ntree=ntree, corr.bias=FALSE)
    RF_model7 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                         .packages = "randomForest") %dopar% randomForest(train_L,WOtrain,ntree=ntree, corr.bias=FALSE)
    RF_model8 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                         .packages = "randomForest") %dopar% randomForest(train_Q,WOtrain,ntree=ntree, corr.bias=FALSE)
    # now we have 4 models, each created using the data in this particular fold
    FM_RF_L[whoop==j,i] <- predict(RF_model3, test_L, type = "response" , predict.all=FALSE)
    FM_RF_Q[whoop==j,i] <- predict(RF_model4, test_Q, type = "response" , predict.all=FALSE)
    WO_RF_L[whoop==j,i] <- predict(RF_model7, test_L, type = "response" , predict.all=FALSE)
    WO_RF_Q[whoop==j,i] <- predict(RF_model8, test_Q, type = "response" , predict.all=FALSE)
  }
}

write.csv(FM_RF_L,"mod3.csv")
write.csv(FM_RF_Q,"mod4.csv")
write.csv(WO_RF_L,"mod7.csv")
write.csv(WO_RF_Q,"mod8.csv")

# save(RF_pred3, RF_pred4, RF_pred7,RF_pred8, file = "RF_pred_v7.rda")

##### Discarded Code ####
# 
# 
# predRF3 <- NULL
# predRF4 <- NULL
# predRF7 <- NULL
# predRF8 <- NULL
# for (i in 1:100)
# {
#   predRF3[[i]] <- rfcv(df1,yFM,ntree=10000, cv.fold=4)
#   predRF4[[i]] <- rfcv(df3,yFM,ntree=10000, cv.fold=4)
#   predRF7[[i]] <- rfcv(df1,yWO,ntree=10000, cv.fold=4)
#   predRF8[[i]] <- rfcv(df3,yWO,ntree=10000, cv.fold=4)
#   print(i)
# }
# for (i in 1:100)
# {
#   # predRF1[[i]] <- cv_RF(feat=df1, out=yFM)
#   # predRF2[[i]] <- cv_RF(feat=df1, out=yPartFM)
#   predRF3[[i]] <- cv_RF(feat=df1, out=yFM)
#   predRF4[[i]] <- cv_RF(feat=df3, out=yFM)
#   # predRF5[[i]] <- cv_RF(feat=df2, out=yPartFM)
#   # predRF6[[i]] <- cv_RF(feat=df2, out=yWO)
#   predRF7[[i]] <- cv_RF(feat=df1, out=yWO)
#   predRF8[[i]] <- cv_RF(feat=df3, out=yWO)
# }
# 
# 
# 
# # RF_pred1 <- cleanup_RF(predRF1)
# # RF_pred2 <- cleanup_RF(predRF2)
# RF_pred3 <- cleanup_RF(predRF3)
# RF_pred4 <- cleanup_RF(predRF4)
# # RF_pred5 <- cleanup_RF(predRF5)
# # RF_pred6 <- cleanup_RF(predRF6)
# RF_pred7 <- cleanup_RF(predRF7)
# RF_pred8 <- cleanup_RF(predRF8)
# predRF3[[i]][whoop==j] <- predict(RF_model3, test_L, type = "response" , predict.all=FALSE)
# predRF4[[i]][whoop==j] <- predict(RF_model4, test_Q, type = "response" , predict.all=FALSE)
# predRF7[[i]][whoop==j] <- predict(RF_model7, test_L, type = "response" , predict.all=FALSE)
# predRF8[[i]][whoop==j] <- predict(RF_model8, test_Q, type = "response" , predict.all=FALSE)
# print("I filled all predictions")
# FM_RF_L[,i] <- predRF3[[i]]
# FM_RF_Q[,i] <- predRF4[[i]]
# WO_RF_L[,i] <- predRF7[[i]]
# WO_RF_Q[,i] <- predRF8[[i]]
