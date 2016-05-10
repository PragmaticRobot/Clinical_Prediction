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

##
Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","lasso_pred.mat")
LASSO_res <- readMat(Dpath, maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","RF_pred.mat")
RF_res <- readMat(Dpath, maxLength=NULL, fixNames=TRUE, Verbose=FALSE)

## LASSO
FM_lin_lasso <-   LASSO_res$LASSO.pred1[[1]]
pFM_lin_lasso <-  LASSO_res$LASSO.pred2[[1]]
WO_lin_lasso  <-  LASSO_res$LASSO.pred3[[1]]
FM_quad_lasso <-  LASSO_res$LASSO.pred4[[1]]
pFM_quad_lasso <- LASSO_res$LASSO.pred5[[1]]
WO_quad_lasso <-  LASSO_res$LASSO.pred6[[1]]

## Random Forests
FM_lin_RF <-   RF_res$RF.pred1[[1]]
pFM_lin_RF <-  RF_res$RF.pred2[[1]]
WO_lin_RF <-   RF_res$RF.pred3[[1]]
FM_quad_RF <-  RF_res$RF.pred4[[1]]
pFM_quad_RF <- RF_res$RF.pred5[[1]]
WO_quad_RF <-  RF_res$RF.pred6[[1]]

## outcomes
yFM <- LASSO_res$yFM
ypFM <- LASSO_res$yPartFM
yWO <- LASSO_res$yWO

## fitting the linear models
fm_lin_lasso_mod <- lm(yFM~rowMeans(FM_lin_lasso))
fm_quad_lasso_mod <- lm(yFM~rowMeans(FM_quad_lasso))
fm_lin_rf_mod <- lm(yFM~rowMeans(FM_lin_RF))
fm_quad_rf_mod <- lm(yFM~rowMeans(FM_quad_RF))

pfm_lin_lasso_mod <- lm(ypFM~rowMeans(pFM_lin_lasso))
pfm_quad_lasso_mod <- lm(ypFM~rowMeans(pFM_quad_lasso))
pfm_lin_rf_mod <- lm(ypFM~rowMeans(pFM_lin_RF))
pfm_quad_rf_mod <- lm(ypFM~rowMeans(pFM_quad_RF))

wo_lin_lasso_mod <- lm(yWO~rowMeans(WO_lin_lasso))
wo_quad_lasso_mod <- lm(yWO~rowMeans(WO_quad_lasso))
wo_lin_rf_mod <- lm(yWO~rowMeans(WO_lin_RF))
wo_quad_rf_mod <- lm(yWO~rowMeans(WO_quad_RF))

## plotting I: Fugl-Meyer, make plot pretty
plot(rowMeans(FM_lin_lasso),yFM, ylim = c(-10, 20), xlim = c(-10, 20), col='red',
     xlab="Predicted Change in Fugl-Meyer Score",ylab="Change in Fugl-Meyer Score",
     cex.lab=0.8,pch=21, cex=0.8)
abline(fm_lin_lasso_mod, col='red')
blah <- summary(fm_lin_lasso_mod)
r1 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(FM_quad_lasso),yFM, col = 'blue',pch=22, cex=0.8)
abline(fm_quad_lasso_mod, col='blue')
blah <- summary(fm_quad_lasso_mod)
r2 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(FM_lin_RF),yFM, col = 'green',pch=24, cex=0.8)
abline(fm_lin_rf_mod, col='green')
blah <- summary(fm_lin_rf_mod)
r3 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(FM_quad_RF),yFM, col = 'orange',pch=25, cex=0.8)
abline(fm_quad_rf_mod, col='orange')
blah <- summary(fm_quad_rf_mod)
r4 <- round(blah$adj.r.squared, digits = 2)
abline(0,1, col='black',pch=16,cex=0.6)
legend("topleft",c(paste('Linear LASSO R^2 = ',as.character(r1)),
                   paste('Quadratic LASSO R^2 = ',as.character(r2)),
                   paste('Linear Random Forests R^2 = ',as.character(r3)),
                   paste('Quadratic Random Forests R^2 = ',as.character(r4))),
       cex=0.8, col=c('red','blue','green','orange'),
       pch=c(21,22,24,25),bty="n")


## plotting II: Arm-only Fugl-Meyer, make plot pretty
plot(rowMeans(pFM_lin_lasso),ypFM, ylim = c(-7, 12), xlim = c(-7, 12), col='red',
     xlab="Predicted Change in Arm Fugl-Meyer Score",ylab="Change in Arm-Only Fugl-Meyer Score",
     cex.lab=0.8,pch=21, cex=0.8)
abline(pfm_lin_lasso_mod, col='red')
blah <- summary(pfm_lin_lasso_mod)
r5 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(pFM_quad_lasso),ypFM, col = 'blue',pch=22, cex=0.8)
abline(pfm_quad_lasso_mod, col='blue')
blah <- summary(pfm_quad_lasso_mod)
r6 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(pFM_lin_RF),ypFM, col = 'green',pch=24, cex=0.8)
abline(pfm_lin_rf_mod, col='green')
blah <- summary(pfm_lin_rf_mod)
r7 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(pFM_quad_RF),ypFM, col = 'orange',pch=25, cex=0.8)
abline(pfm_quad_rf_mod, col='orange')
blah <- summary(pfm_quad_rf_mod)
r8 <- round(blah$adj.r.squared, digits = 2)
abline(0,1, col='black',pch=16,cex=0.6)
legend("topleft",c(paste('Linear LASSO R^2 = ',as.character(r5)),
                   paste('Quadratic LASSO R^2 = ',as.character(r6)),
                   paste('Linear Random Forests R^2 = ',as.character(r7)),
                   paste('Quadratic Random Forests R^2 = ',as.character(r8))),
       cex=0.8, col=c('red','blue','green','orange'),
       pch=c(21,22,24,25),bty="n")


## plotting III: Wolf, make plot pretty
plot(rowMeans(WO_lin_lasso),yWO, ylim = c(-28, 14), xlim = c(-28, 14), col='red',
     xlab="Predicted Change in Wolf Score",ylab="Change in Wolf Motor Function Score",
     cex.lab=0.8,pch=21, cex=0.8)
abline(wo_lin_lasso_mod, col='red')
blah <- summary(wo_lin_lasso_mod)
r9 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(WO_quad_lasso),yWO, col = 'blue',pch=22, cex=0.8)
abline(wo_quad_lasso_mod, col='blue')
blah <- summary(wo_quad_lasso_mod)
r10 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(WO_lin_RF),yWO, col = 'green',pch=24, cex=0.8)
abline(wo_lin_rf_mod, col='green')
blah <- summary(wo_lin_rf_mod)
r11 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(WO_quad_RF),yWO, col = 'orange',pch=25, cex=0.8)
abline(wo_quad_rf_mod, col='orange')
blah <- summary(wo_quad_rf_mod)
r12 <- round(blah$adj.r.squared, digits = 2)
abline(0,1, col='black',pch=16,cex=0.6)
legend("topleft",c(paste('Linear LASSO R^2 = ',as.character(r9)),
                   paste('Quadratic LASSO R^2 = ',as.character(r10)),
                   paste('Linear Random Forests R^2 = ',as.character(r11)),
                   paste('Quadratic Random Forests R^2 = ',as.character(r12))),
       cex=0.8, col=c('red','blue','green','orange'),
       pch=c(21,22,24,25),bty="n")

