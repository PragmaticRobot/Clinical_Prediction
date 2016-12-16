# This version adds bootstrapping to LASSO, so the model predictions
# can be averaged and mean and std can be calculated for these
# predictions, then we plot the fits and calculate R^2 (coef of determination)

rm(list = ls()) # clear the environment

#################################
########## Libraries ############
#################################
# options(error=recover)
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,e1071,
               psych, foreach, randomForest, doParallel, inTrees,tableone,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)
source("functions.R")
source("create_df1_df2_no_unknown.R")

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations


##
load('Basics.rda')
# Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","lasso_pred.mat")
# LASSO_res <- readMat("lasso_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
load("2016dec12_LS_nReps1_nFolds26.rda")
# Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","RF_pred.mat")
# RF_res <- readMat("RF_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
load("2016dec12_RF_nReps1_nFolds26.rda")
# 
# ## LASSO
# # FM_lin_lasso <-   LASSO_res$LASSO.pred1[[1]]
# FM_lin_lasso <-   LASSO_pred1[[1]]
# # pFM_lin_lasso <-  LASSO_res$LASSO.pred2[[1]]
# # WO_lin_lasso  <-  LASSO_res$LASSO.pred3[[1]]
# WO_lin_lasso  <-  LASSO_pred5[[1]]
# # FM_quad_lasso <-  LASSO_res$LASSO.pred4[[1]]
# FM_quad_lasso <-  LASSO_pred2[[1]]
# # pFM_quad_lasso <- LASSO_res$LASSO.pred5[[1]]
# # WO_quad_lasso <-  LASSO_res$LASSO.pred6[[1]]
# WO_quad_lasso <-  LASSO_pred6[[1]]
# 
# ## Random Forests
# # FM_lin_RF <-   RF_res$RF.pred1[[1]]
# FM_lin_RF <-   RF_pred3[[1]]
# # pFM_lin_RF <-  RF_res$RF.pred2[[1]]
# # WO_lin_RF <-   RF_res$RF.pred3[[1]]
# WO_lin_RF <-   RF_pred7[[1]]
# # FM_quad_RF <-  RF_res$RF.pred4[[1]]
# FM_quad_RF <-  RF_pred4[[1]]
# # pFM_quad_RF <- RF_res$RF.pred5[[1]]
# # WO_quad_RF <-  RF_res$RF.pred6[[1]]
# WO_quad_RF <-  RF_pred8[[1]]

## outcomes
# yFM <- LASSO_res$yFM
# ypFM <- LASSO_res$yPartFM
# yWO <- LASSO_res$yWO

## fitting the linear models
fm_lin_lasso_mod <- lm(yFM~rowMeans(FM_LS_L))
fm_quad_lasso_mod <- lm(yFM~rowMeans(FM_LS_Q))
fm_lin_rf_mod <- lm(yFM~rowMeans(FM_RF_L))
fm_quad_rf_mod <- lm(yFM~rowMeans(FM_RF_Q))
# fm_noextrap_lasso_mod <- lm(rowMeans(FM_NoExtrap_lasso)~yFM)
# fm_noextrap_rf_mod <- lm(rowMeans(FM_NoExtrap_RF)~yFM)

# pfm_lin_lasso_mod <- lm(rowMeans(pFM_lin_lasso)~ypFM)
# pfm_quad_lasso_mod <- lm(rowMeans(pFM_quad_lasso)~ypFM)
# pfm_lin_rf_mod <- lm(rowMeans(pFM_lin_RF)~ypFM)
# pfm_quad_rf_mod <- lm(rowMeans(pFM_quad_RF)~ypFM)

wo_lin_lasso_mod <- lm(yWO~rowMeans(WO_LS_L))
wo_quad_lasso_mod <- lm(yWO~rowMeans(WO_LS_Q))
wo_lin_rf_mod <- lm(yWO~rowMeans(WO_RF_L))
wo_quad_rf_mod <- lm(yWO~rowMeans(WO_RF_Q))
# wo_noextrap_lasso_mod <- lm(rowMeans(WO_NoExtrap_lasso)~yWO)
# wo_noextrap_rf_mod <- lm(rowMeans(WO_NoExtrap_RF)~yWO)

# put all the inputs (rowMeans) in a matrix for passing to the plotting function

InsFM <- cbind(rowMeans(FM_LS_L),rowMeans(FM_LS_Q),
                 rowMeans(FM_RF_L), rowMeans(FM_RF_Q))
               # rowMeans(FM_NoExtrap_lasso), rowMeans(FM_NoExtrap_RF)) # this is for UEFM

# InsPFM <- cbind(rowMeans(pFM_lin_lasso),rowMeans(pFM_quad_lasso),
#                rowMeans(pFM_lin_RF), rowMeans(pFM_quad_RF)) # this is for Arm FM

InsWO <- cbind(rowMeans(WO_LS_L),rowMeans(WO_LS_Q),
                 rowMeans(WO_RF_L), rowMeans(WO_RF_Q))
               # rowMeans(WO_NoExtrap_lasso), rowMeans(WO_NoExtrap_RF)) # this is for Wolf

modListFM <- list(fm_lin_lasso_mod,fm_quad_lasso_mod,
                  fm_lin_rf_mod,fm_quad_rf_mod)
                  # fm_noextrap_lasso_mod,fm_noextrap_rf_mod)
# modListPFM <- list(pfm_lin_lasso_mod,pfm_quad_lasso_mod,
                  # pfm_lin_rf_mod,pfm_quad_rf_mod)
modListWO <- list(wo_lin_lasso_mod,wo_quad_lasso_mod,
                  wo_lin_rf_mod,wo_quad_rf_mod)
                  # wo_noextrap_lasso_mod, wo_noextrap_rf_mod)

# last input to following function:
# 1: Use RMSE
# 2: Use Slope (correlation)
# 3: Use Coef of Determination (R^2)
plot_LSRF_fits(InsFM,modListFM,relative=0,yFM,'Fugl-Meyer',1)
# plot_LSRF_fits(InsPFM,modListPFM,relative=0,ypFM,'Arm-Only Fugl-Meyer')
plot_LSRF_fits(InsWO,modListWO,relative=0,yWO,'Wolf',1)

## Create and use a DoHists function here

skFM <- DoHists(InsFM,yFM,'Fugl-Meyer')
# skpFM <- DoHists(InsPFM,ypFM,'Arm-Only Fugl-Meyer')
skWO <- DoHists(InsWO,yWO,'Wolf')

skFM2 <- DoLines(InsFM,yFM,'Fugl-Meyer')
skWO2 <- DoLines(InsWO,yWO,'Wolf')
# skpFM <- DoHists(InsPFM,ypFM,'Arm-Only Fugl-Meyer')

# write.csv(FM_lin_lasso,"mod1.csv")
# # write.csv(FM_NoExtrap_lasso,"mod2.csv")
# write.csv(FM_quad_lasso,"mod2.csv")
# write.csv(FM_lin_RF,"mod3.csv")
# # write.csv(FM_NoExtrap_RF,"mod5.csv")
# write.csv(FM_quad_RF,"mod4.csv")
# 
# write.csv(WO_lin_lasso,"mod5.csv")
# # write.csv(WO_NoExtrap_lasso,"mod8.csv")
# write.csv(WO_quad_lasso,"mod6.csv")
# write.csv(WO_lin_RF,"mod7.csv")
# # write.csv(WO_NoExtrap_RF,"mod11.csv")
# write.csv(WO_quad_RF,"mod8.csv")
# 
# write.csv(yFM,"yFM.csv")
# write.csv(yWO,"yWO.csv")
