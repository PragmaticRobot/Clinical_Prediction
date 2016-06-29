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
registerDoParallel(cores=4) # 4 cores to do the simulations

load("colnames.rda")

################ The BIG Model (all 55 features) ############

predRF1 <- NULL
predRF2 <- NULL

for (i in 1:100)
{
  predRF1[[i]] <- cv_RF(feat=Features, out=yFM)
  predRF2[[i]] <- cv_RF(feat=Features, out=yWO)
}

RF_pred1 <- cleanup_RF(predRF1)
RF_pred2 <- cleanup_RF(predRF2)

RF_1 <- data.frame(matrix(vector(),55,100), row.names = gverb)
RF_2 <- data.frame(matrix(vector(),55,100), row.names = gverb)

for (j in 1:100) # loop over cross-validations
{
  for (i in 1:55)
  {
    RF_1[RF_pred1[[2]][i,j],j] <- RF_pred1[[3]][i,j]
    RF_2[RF_pred2[[2]][i,j],j] <- RF_pred2[[3]][i,j]
  }
}

RF1_Ind <- RF_1[order(-RF_1$X1), , drop=FALSE]
RF2_Ind <- RF_2[order(-RF_2$X1), , drop=FALSE]
Ind_RF_FM <- order(-RF_1$X1)
Ind_RF_WO <- order(-RF_2$X1)


save(RF_1,RF_2,file = "RF_preds.rda")


## Find original RMSE/R^2/Slope
pXX <- as.vector(predict_Full)
# RMSE
full_rmse <- round(sqrt(mean((yFM-pXX)^2)), digits = 2)
mod_full <- lm(pXX~yFM)
blah <- summary(mod_full)
full_r2 <- round(blah$adj.r.squared, digits = 2)
blah <- summary(mod_full,correlation=TRUE)
blah <- round(blah$correlation, digits = 2)
full_sl <- blah[2]

############ Initialize the hit ##########
hit_rmse <- matrix(0, 55, 1)
hit_r2 <- matrix(0, 55, 1)
hit_sl <- matrix(0, 55, 1)
## reorder df1 by variable importance returned according to Index (line 29)
dfFM_VI <- df1[, order(Ind_RF_FM),drop = FALSE]
dfWO_VI <- df1[, order(Ind_RF_WO),drop = FALSE]

# df_VI now contains reordered such that the first column is the most
# important feature, and the last column is the least important
hit_rmse[55] <- full_rmse
hit_r2[55] <- full_r2
hit_sl[55] <- full_sl
# now the last element of the hit vectors contains values for the full model
# finally we can start omitting variables (starting with the least important)
# and see the effect on the hit list

save(dfFM_VI,dfWO_VI,hit_sl,hit_r2,hit_rmse,file = "df_VI_hitlist.rda")