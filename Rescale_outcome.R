# This script will rescale the outcome (change in clinical score) to try to account for
# non-normality as per Konrad's suggestion. This reconstructs the models, then plots the new
# fits and histograms.

# Two transformations of the outcome are suggested: relative change in score
# and log(change in score)
# just because, will not do Arm-only FM

#################################
########## Libraries ############
#################################
# options(error=recover)
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,
               psych, foreach, randomForest, doParallel, inTrees,tableone,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)
source("functions.R")
source("create_df1_df2.R")

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations

## Log change in score doesn't work, some changes are negative, so we get NaNs
# lgFM <- log(yFM)
# lgWO <- log(yWO)

rvFM <- yFM/df1$init_FM
rvWO <- yWO/df1$initWO

predLASSO_rv1<-NULL
predLASSO_rv2<-NULL
predLASSO_rv3<-NULL
predLASSO_rv4<-NULL

## Create a progress bar to monitor this shit:
pb1 <- tkProgressBar(title = "LASSO models", min = 0, max = 100, width = 400)
for (i in 1:100)
{
  predLASSO_rv1[[i]] <- cv_mod(feat=df1,out=rvFM,stand=1,isWO = 0)
  predLASSO_rv2[[i]] <- cv_mod(feat=df1,out=rvWO,stand=1,isWO = 1)
  predLASSO_rv3[[i]] <- cv_mod(feat=df2,out=rvFM,stand=1,isWO = 0)
  predLASSO_rv4[[i]] <- cv_mod(feat=df2,out=rvWO,stand=1,isWO = 1)
  setTkProgressBar(pb1,i,label=paste(round(i/100*100,0),"% done"))
}
close(pb1)

LASSO_pred_rv1 <- cleanup_LASSO(predLASSO_rv1)
LASSO_pred_rv2 <- cleanup_LASSO(predLASSO_rv2)
LASSO_pred_rv3 <- cleanup_LASSO(predLASSO_rv3)
LASSO_pred_rv4 <- cleanup_LASSO(predLASSO_rv4)

predRF1 <- NULL
predRF2 <- NULL
predRF3 <- NULL
predRF4 <- NULL

pb2 <- tkProgressBar(title = "Random Forest models", min = 0, max = 100, width = 400)

for (i in 1:100)
{
  predRF1[[i]] <- cv_RF(feat=df1, out=rvFM)
  predRF2[[i]] <- cv_RF(feat=df1, out=rvWO)
  predRF3[[i]] <- cv_RF(feat=df2, out=rvFM)
  predRF4[[i]] <- cv_RF(feat=df2, out=rvWO)
  setTkProgressBar(pb2,i,label=paste(round(i/100*100,0),"% done"))
}
close(pb2)

RF_pred_rv1 <- cleanup_RF(predRF1)
RF_pred_rv2 <- cleanup_RF(predRF2)
RF_pred_rv3 <- cleanup_RF(predRF3)
RF_pred_rv4 <- cleanup_RF(predRF4)

## LASSO
rvFM_lin_lasso <- LASSO_pred_rv1[[1]]
rvWO_lin_lasso <- LASSO_pred_rv2[[1]]
rvFM_quad_lasso <- LASSO_pred_rv3[[1]]
rvWO_quad_lasso <- LASSO_pred_rv4[[1]]

## Random Fores1ts
rvFM_lin_RF <- RF_pred_rv1[[1]]
rvWO_lin_RF <- RF_pred_rv2[[1]]
rvFM_quad_RF <- RF_pred_rv3[[1]]
rvWO_quad_RF <- RF_pred_rv4[[1]]

save(rvFM_lin_lasso, rvFM_quad_lasso, rvFM_lin_RF, rvFM_quad_RF, rvWO_lin_lasso,
     rvWO_quad_lasso, rvWO_lin_RF, rvWO_quad_RF, file = "rv_cv_models.rda")

save(RF_pred_rv1,RF_pred_rv2,RF_pred_rv3,RF_pred_rv4,LASSO_pred_rv1,LASSO_pred_rv2,
     LASSO_pred_rv3,LASSO_pred_rv4,file = "rv_cv_all_LSRF.rda")

save(rvFM,rvWO, file = "rescaled_out.rda")
