## This script does backward selection on features picked by random forests
# and LASSO (later) and evaluates the effects on three outcomes:
# 1- RMSE
# 2- Coef of Determination
# 3- Slope

par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,
               psych, foreach, randomForest, doParallel, inTrees,tableone,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)
source("functions.R")
source("create_df1_df2.R")

## set up the parallel pool
registerDoParallel(cores=8) # 8 cores to do the simulations
RF_res <- readMat("RF_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)

load("colnames.rda")

################ The BIG Model (all 55 features) ############

RF_1 <- data.frame(matrix(vector(),55,100), row.names = gverb)
RF_2 <- data.frame(matrix(vector(),55,100), row.names = gverb)

for (j in 1:100) # loop over cross-validations
{
  for (i in 1:55)
  {
    RF_1[RF_res$RF.pred1[[2]][i,j],j] <- RF_res$RF.pred1[[3]][i,j]
    RF_2[RF_res$RF.pred3[[2]][i,j],j] <- RF_res$RF.pred3[[3]][i,j]
  }
}
RF_FM_pred <- RF_res$RF.pred1[[1]]
RF_WO_pred <- RF_res$RF.pred3[[1]]
RF1_Ind <- RF_1[order(-RF_1$X1), , drop=FALSE]
RF2_Ind <- RF_2[order(-RF_2$X1), , drop=FALSE]
Ind_RF_FM <- order(-RF_1$X1)
Ind_RF_WO <- order(-RF_2$X1)

save(RF_1,RF_2,RF1_Ind,RF2_Ind,file = "RF_preds.rda")

############ Initialize the hit ##########
## Find original RMSE/R^2/Slope
pXX <- rowMeans(RF_FM_pred)
pyy <- rowMeans(RF_WO_pred)
hit_rmse <- matrix(0, 55, 1)
hit_r2 <- matrix(0, 55, 1)
hit_sl <- matrix(0, 55, 1)
hitlist_FM <- data.frame(hit_rmse,hit_r2,hit_sl)
hitlist_WO <- data.frame(hit_rmse,hit_r2,hit_sl)
colnames(hitlist_FM) <- c("RMSE","Coef Det.","Slope")
colnames(hitlist_WO) <- c("RMSE","Coef Det.","Slope")


hitlist_FM$RMSE[55] <- round(sqrt(mean((yFM-pXX)^2)), digits = 2)
mod_full <- lm(pXX~yFM)
blah <- summary(mod_full)
hitlist_FM$`Coef Det.`[55] <- round(blah$adj.r.squared, digits = 2)
blah <- summary(mod_full,correlation=TRUE)
blah <- round(blah$correlation, digits = 2)
hitlist_FM$Slope[55] <- blah[2]

hitlist_WO$RMSE[55] <- round(sqrt(mean((yWO-pyy)^2)), digits = 2)
mod_full <- lm(pyy~yWO)
blah <- summary(mod_full)
hitlist_WO$`Coef Det.`[55] <- round(blah$adj.r.squared, digits = 2)
blah <- summary(mod_full,correlation=TRUE)
blah <- round(blah$correlation, digits = 2)
hitlist_WO$Slope[55] <- blah[2]

## reorder df1 by variable importance returned according to Index (line 29)
dfFM_VI <- df1[, Ind_RF_FM,drop = FALSE]
dfWO_VI <- df1[, Ind_RF_WO,drop = FALSE]

# df_VI now contains reordered such that the first column is the most
# important feature, and the last column is the least important

# now the last element of the hit vectors contains values for the full model
# finally we can start omitting variables (starting with the least important)
# and see the effect on the hit list

save(dfFM_VI,dfWO_VI,hit_sl,hit_r2,hit_rmse,file = "df_VI_hitlist.rda")

############## Start deleting least important features ###########
gverbInd <- gverb[Ind_RF_FM]
dfFM <- dfFM_VI
dfWO <- dfWO_VI
for (i in 1:53){
  ind <- 55-i
  dfFM <- dfFM[,-c(ind+1)]
  dfWO <- dfWO[,-c(ind+1)]
  gverbInd <- gverbInd[-c(ind+1)]
  predRF1 <- NULL
  predRF2 <- NULL
  predRF1 <- rf_just_pred(feat=dfFM, out=yFM)
  predRF2 <- rf_just_pred(feat=dfWO, out=yWO)
  # RF_pred1 <- cleanup_RF(predRF1)
  # RF_pred2 <- cleanup_RF(predRF2)
  
  # RF_FM_pred <- RF_pred1[[1]]
  # RF_WO_pred <- RF_pred2[[1]]
  RF_FM_pred <- as.vector(predRF1)
  RF_WO_pred <- as.vector(predRF2)
  # pXX <- rowMeans(RF_FM_pred)
  # pyy <- rowMeans(RF_WO_pred)
  pXX <- RF_FM_pred
  pyy <- RF_WO_pred
 
  hitlist_FM$RMSE[ind] <- round(sqrt(mean((yFM-pXX)^2)), digits = 2)
  mod_full <- lm(pXX~yFM)
  blah <- summary(mod_full)
  hitlist_FM$`Coef Det.`[ind] <- round(blah$adj.r.squared, digits = 2)
  blah <- summary(mod_full,correlation=TRUE)
  blah <- round(blah$correlation, digits = 2)
  hitlist_FM$Slope[ind] <- blah[2]
  
  hitlist_WO$RMSE[ind] <- round(sqrt(mean((yWO-pyy)^2)), digits = 2)
  mod_full <- lm(pyy~yWO)
  blah <- summary(mod_full)
  hitlist_WO$`Coef Det.`[ind] <- round(blah$adj.r.squared, digits = 2)
  blah <- summary(mod_full,correlation=TRUE)
  blah <- round(blah$correlation, digits = 2)
  hitlist_WO$Slope[ind] <- blah[2]
}
hitlist_WO <- hitlist_WO[-c(1),]
hitlist_FM <- hitlist_FM[-c(1),]

save(hitlist_FM,hitlist_WO, file = "hitlist.rda")

plot(df1$init_FM,yFM, xlab="Initial UEFM",ylab="Change in UEFM")
abline(lm(yFM~df1$init_FM))
plot(df1$MeanMaxSpd,yFM, xlab="Mean Max Speed",ylab="Change In UEFM")
abline(lm(yFM~df1$MeanMaxSpd))
