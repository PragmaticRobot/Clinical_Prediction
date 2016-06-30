# This version adds bootstrapping to LASSO, so the model predictions
# can be averaged and mean and std can be calculated for these
# predictions, then we plot the fits and calculate R^2 (coef of determination)

# rm(list = ls()) # clear the environment

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
source("create_df1_df2.R")

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations


##
Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","lasso_pred.mat")
LASSO_res <- readMat(Dpath, maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","RF_pred.mat")
RF_res <- readMat(Dpath, maxLength=NULL, fixNames=TRUE, Verbose=FALSE)

## LASSO
FM_lin_lasso <-   LASSO_res$LASSO.pred1[[1]]
# pFM_lin_lasso <-  LASSO_res$LASSO.pred2[[1]]
WO_lin_lasso  <-  LASSO_res$LASSO.pred3[[1]]
FM_quad_lasso <-  LASSO_res$LASSO.pred4[[1]]
# pFM_quad_lasso <- LASSO_res$LASSO.pred5[[1]]
WO_quad_lasso <-  LASSO_res$LASSO.pred6[[1]]

## Random Forests
FM_lin_RF <-   RF_res$RF.pred1[[1]]
# pFM_lin_RF <-  RF_res$RF.pred2[[1]]
WO_lin_RF <-   RF_res$RF.pred3[[1]]
FM_quad_RF <-  RF_res$RF.pred4[[1]]
# pFM_quad_RF <- RF_res$RF.pred5[[1]]
WO_quad_RF <-  RF_res$RF.pred6[[1]]

## outcomes
yFM <- LASSO_res$yFM
# ypFM <- LASSO_res$yPartFM
yWO <- LASSO_res$yWO

## fitting the linear models
fm_lin_lasso_mod <- lm(rowMeans(FM_lin_lasso)~yFM)
fm_quad_lasso_mod <- lm(rowMeans(FM_quad_lasso)~yFM)
fm_lin_rf_mod <- lm(rowMeans(FM_lin_RF)~yFM)
fm_quad_rf_mod <- lm(rowMeans(FM_quad_RF)~yFM)

# pfm_lin_lasso_mod <- lm(rowMeans(pFM_lin_lasso)~ypFM)
# pfm_quad_lasso_mod <- lm(rowMeans(pFM_quad_lasso)~ypFM)
# pfm_lin_rf_mod <- lm(rowMeans(pFM_lin_RF)~ypFM)
# pfm_quad_rf_mod <- lm(rowMeans(pFM_quad_RF)~ypFM)

wo_lin_lasso_mod <- lm(rowMeans(WO_lin_lasso)~yWO)
wo_quad_lasso_mod <- lm(rowMeans(WO_quad_lasso)~yWO)
wo_lin_rf_mod <- lm(rowMeans(WO_lin_RF)~yWO)
wo_quad_rf_mod <- lm(rowMeans(WO_quad_RF)~yWO)


# put all the inputs (rowMeans) in a matrix for passing to the plotting function

InsFM <- cbind(rowMeans(FM_lin_lasso),rowMeans(FM_quad_lasso),
                 rowMeans(FM_lin_RF), rowMeans(FM_quad_RF)) # this is for UEFM

# InsPFM <- cbind(rowMeans(pFM_lin_lasso),rowMeans(pFM_quad_lasso),
#                rowMeans(pFM_lin_RF), rowMeans(pFM_quad_RF)) # this is for Arm FM

InsWO <- cbind(rowMeans(WO_lin_lasso),rowMeans(WO_quad_lasso),
                 rowMeans(WO_lin_RF), rowMeans(WO_quad_RF)) # this is for Wolf

modListFM <- list(fm_lin_lasso_mod,fm_quad_lasso_mod,
                  fm_lin_rf_mod,fm_quad_rf_mod)
# modListPFM <- list(pfm_lin_lasso_mod,pfm_quad_lasso_mod,
                  # pfm_lin_rf_mod,pfm_quad_rf_mod)
modListWO <- list(wo_lin_lasso_mod,wo_quad_lasso_mod,
                  wo_lin_rf_mod,wo_quad_rf_mod)

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