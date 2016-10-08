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

# code:
# FM: Fugl-Myer
# WO: Wolf Motor Function
# LS: LASSO
# RF: Random FOrests
# L: Linear
# Q: Quadratic, no extrapolation

# Set the cv parameters
nreps <- 1
nfolds <- 26
nSubj <- 26

# to save the predictions
FM_LS_L <- matrix(0,nSubj,nreps)
FM_LS_Q <- matrix(0,nSubj,nreps)
WO_LS_L <- matrix(0,nSubj,nreps)
WO_LS_Q <- matrix(0,nSubj,nreps)
# to save the winning models. a0 is intercept, beta is coefficients
FM_LS_L_a0 <- matrix(0,nreps,nfolds)
FM_LS_L_beta <- rep( list(list()), nreps )

FM_LS_Q_a0 <- matrix(0,nreps,nfolds)
FM_LS_Q_beta <- rep( list(list()), nreps )

WO_LS_L_a0 <- matrix(0,nreps,nfolds)
WO_LS_L_beta <- rep( list(list()), nreps )

WO_LS_Q_a0 <- matrix(0,nreps,nfolds)
WO_LS_Q_beta <- rep( list(list()), nreps )

for (i in 1:nreps)
{
  print(i)
  whoop <- createFolds(yFM, k=nfolds, list=FALSE)
  for (j in 1:nfolds) # 4-fold cross-validation
  {
    # first let's form the train and test sets
    train_L <- df1[whoop!=j,]; test_L <- df1[whoop==j,]; FMtrain <- yFM[whoop!=j]; WOtrain <- yWO[whoop!=j];
    train_Q <- df3[whoop!=j,]; test_Q <- df3[whoop==j,]; FMtest <- yFM[whoop==j]; WOtest <- yWO[whoop==j];
    
    trash <- glmnet(data.matrix(train_L),FMtrain, family = c("gaussian"),alpha = 1, nlambda = 1000,standardize = TRUE)
    trash_mod_ind <-  min(which(trash$dev.ratio > 0.95))
    FM_LS_L[whoop==j,i] <- predict(trash,s=trash$lambda[trash_mod_ind],newx = data.matrix(test_L))
    FM_LS_L_a0[i,j] <- trash$a0[[trash_mod_ind]]
    FM_LS_L_beta[[i]][[j]] <- trash$beta[,trash_mod_ind]
    
    trash <- glmnet(data.matrix(train_Q),FMtrain, family = c("gaussian"),alpha = 1, nlambda = 1000,standardize = TRUE)
    trash_mod_ind <-  min(which(trash$dev.ratio > 0.95))
    FM_LS_Q[whoop==j,i] <- predict(trash,s=trash$lambda[trash_mod_ind],newx = data.matrix(test_Q))
    FM_LS_Q_a0[i,j] <- trash$a0[[trash_mod_ind]]
    FM_LS_Q_beta[[i]][[j]] <- trash$beta[,trash_mod_ind]
    
    trash <- glmnet(data.matrix(train_L),WOtrain, family = c("gaussian"),alpha = 1, nlambda = 1000,standardize = TRUE)
    trash_mod_ind <-  min(which(trash$dev.ratio > 0.95))
    WO_LS_L[whoop==j,i] <- predict(trash,s=trash$lambda[trash_mod_ind],newx = data.matrix(test_L))
    WO_LS_L_a0[i,j] <- trash$a0[[trash_mod_ind]]
    WO_LS_L_beta[[i]][[j]] <- trash$beta[,trash_mod_ind]
    
    trash <- glmnet(data.matrix(train_Q),WOtrain, family = c("gaussian"),alpha = 1, nlambda = 1000,standardize = TRUE)
    trash_mod_ind <-  min(which(trash$dev.ratio > 0.95))
    WO_LS_Q[whoop==j,i] <- predict(trash,s=trash$lambda[trash_mod_ind],newx = data.matrix(test_Q))
    WO_LS_Q_a0[i,j] <- trash$a0[[trash_mod_ind]]
    WO_LS_Q_beta[[i]][[j]] <- trash$beta[,trash_mod_ind]
    
  }
}

write.csv(FM_LS_L,"mod1_v9.csv")
write.csv(FM_LS_Q,"mod2_v9.csv")
write.csv(WO_LS_L,"mod5_v9.csv")
write.csv(WO_LS_Q,"mod6_v9.csv")

save(FM_LS_L_a0, FM_LS_L_beta, FM_LS_Q_a0,FM_LS_Q_beta,WO_LS_L_a0,WO_LS_L_beta,WO_LS_Q_a0,WO_LS_Q_beta,
     FM_LS_L,FM_LS_Q,WO_LS_L,WO_LS_Q, file = "LS_v9.rda")

# #######################################
# #### 06 - Random Forests Bootstrapping ####
# #######################################
# # 3: linear FM RF
# # 4: quad df3 FM RF
# # 7: linear WO RF
# # 8: quad df3 WO RF
# # Set the cross-validation parameters
# nreps <- 1 # number of cross-validation repeats
# nfolds <- 26 # number of folds in each cv
# nSubj <- 26 # number of observations to predict (number of subjects)
# 
# # save the predictions
# FM_RF_L <- matrix(0,nSubj,nreps)
# FM_RF_Q <- matrix(0,nSubj,nreps)
# WO_RF_L <- matrix(0,nSubj,nreps)
# WO_RF_Q <- matrix(0,nSubj,nreps)
# 
# # save the feature importance ranks
# FM_RF_L_VI <- rep( list(list()), nreps )
# FM_RF_Q_VI <- rep( list(list()), nreps )
# WO_RF_L_VI <- rep( list(list()), nreps )
# WO_RF_Q_VI <- rep( list(list()), nreps )
# 
# for (i in 1:nreps)
# {
#   print(i)
#   whoop <- createFolds(yFM, k=nfolds, list=FALSE)
#   for (j in 1:nfolds) # nfolds-fold cross-validation
#   {
#     # first let's form the train and test sets
#     train_L <- df1[whoop!=j,]; test_L <- df1[whoop==j,]; FMtrain <- yFM[whoop!=j]; WOtrain <- yWO[whoop!=j];
#     train_Q <- df3[whoop!=j,]; test_Q <- df3[whoop==j,]; FMtest <- yFM[whoop==j]; WOtest <- yWO[whoop==j];
#     RF_model3 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
#                          .packages = "randomForest") %dopar% randomForest(train_L,FMtrain,ntree=ntree, corr.bias=FALSE)
#     FM_RF_L_VI[[i]][[j]] <- RF_model3$importance
# 
#     RF_model4 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
#                          .packages = "randomForest") %dopar% randomForest(train_Q,FMtrain,ntree=ntree, corr.bias=FALSE)
#     FM_RF_Q_VI[[i]][[j]] <- RF_model4$importance
#     
#     RF_model7 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
#                          .packages = "randomForest") %dopar% randomForest(train_L,WOtrain,ntree=ntree, corr.bias=FALSE)
#     WO_RF_L_VI[[i]][[j]] <- RF_model7$importance
#     
#     RF_model8 <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
#                          .packages = "randomForest") %dopar% randomForest(train_Q,WOtrain,ntree=ntree, corr.bias=FALSE)
#     WO_RF_Q_VI[[i]][[j]] <- RF_model8$importance
#     
#     # now we have 4 models, each created using the data in this particular fold
#     FM_RF_L[whoop==j,i] <- predict(RF_model3, test_L, type = "response" , predict.all=FALSE)
#     FM_RF_Q[whoop==j,i] <- predict(RF_model4, test_Q, type = "response" , predict.all=FALSE)
#     WO_RF_L[whoop==j,i] <- predict(RF_model7, test_L, type = "response" , predict.all=FALSE)
#     WO_RF_Q[whoop==j,i] <- predict(RF_model8, test_Q, type = "response" , predict.all=FALSE)
#   }
# }
# 
# write.csv(FM_RF_L,"mod3_v9.csv")
# write.csv(FM_RF_Q,"mod4_v9.csv")
# write.csv(WO_RF_L,"mod7_v9.csv")
# write.csv(WO_RF_Q,"mod8_v9.csv")
# 
# save(FM_RF_L,FM_RF_Q,WO_RF_L,WO_RF_Q,FM_RF_L_VI,FM_RF_Q_VI,WO_RF_L_VI,WO_RF_Q_VI, file = "RF_v9.rda")