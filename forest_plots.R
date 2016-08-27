#### Prep ####

# par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,forestplot,
               psych, foreach, randomForest, doParallel, inTrees,tableone,boot,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)
source("functions.R")
source("create_df1_df2.R")

registerDoParallel(cores=4) # 4 cores to do the simulations

load("LS_RF_Results.rda")

# Order of variables 1-8
# FM_Lin_LASSO_sort
# FM_Lin_RF_sort
# FM_Quad_LASSO_sort
# FM_Quad_RF_sort
# Wolf_Lin_LASSO_sort
# Wolf_Lin_RF_sort
# Wolf_Quad_LASSO_sort
# Wolf_Quad_RF_sort

##### Initialize empty vectors to hold data (LowCI, highCI) ####
low1 <- vector()
low2 <- vector()
low3 <- vector()
low4 <- vector()
low5 <- vector()
low6 <- vector()
low7 <- vector()
low8 <- vector()

high1 <- vector()
high2 <- vector()
high3 <- vector()
high4 <- vector()
high5 <- vector()
high6 <- vector()
high7 <- vector()
high8 <- vector()

#### Populate the low/high vectors ####
for (i in 1:length(Booted1)){
  low1[i] <- as.double(Booted1[[i]][1])
  high1[i] <- as.double(Booted1[[i]][2])
}

for (i in 1:length(Booted2)){
  low2[i] <- as.double(Booted2[[i]][1])
  high2[i] <- as.double(Booted2[[i]][2])
}

for (i in 1:length(Booted3)){
  low3[i] <- as.double(Booted3[[i]][1])
  high3[i] <- as.double(Booted3[[i]][2])
}

for (i in 1:length(Booted4)){
  low4[i] <- as.double(Booted4[[i]][1])
  high4[i] <- as.double(Booted4[[i]][2])
}

for (i in 1:length(Booted5)){
  low5[i] <- as.double(Booted5[[i]][1])
  high5[i] <- as.double(Booted5[[i]][2])
}

for (i in 1:length(Booted6)){
  low6[i] <- as.double(Booted6[[i]][1])
  high6[i] <- as.double(Booted6[[i]][2])
}

for (i in 1:length(Booted7)){
  low7[i] <- as.double(Booted7[[i]][1])
  high7[i] <- as.double(Booted7[[i]][2])
}

for (i in 1:length(Booted8)){
  low8[i] <- as.double(Booted8[[i]][1])
  high8[i] <- as.double(Booted8[[i]][2])
}

#### Extract variable names from LS_RF_Results for plotting ####
rnam1 <- rownames(FM_Lin_LASSO_sort)
rnam2 <- rownames(FM_Lin_RF_sort)
rnam3 <- rownames(FM_Quad_LASSO_sort)
rnam4 <- rownames(FM_Quad_RF_sort)
rnam5 <- rownames(Wolf_Lin_LASSO_sort)
rnam6 <- rownames(Wolf_Lin_RF_sort)
rnam7 <- rownames(Wolf_Quad_LASSO_sort)
rnam8 <- rownames(Wolf_Quad_RF_sort)

#### Get means and medians of feature sets ####

mns1 <- apply(as.matrix(df1), 1, mean) 
mns2 <- apply(as.matrix(df2), 1, mean)
meds1 <- apply(as.matrix(df1), 1, median)
meds2 <- apply(as.matrix(df2), 1, median)

#### Notes at this point: ####
# df1 contains original feature set
# df2 contains quadratic feature set
# we have means, medians, lowCIs, highCIs
# we use the means for scaling purposes
# forest plot the medians and CIs

#### Forest plot ####

## Interwebz code:
# library(forestplot)
# # Cochrane data from the 'rmeta'-package
# cochrane_from_rmeta <- 
#   structure(list(
#     mean  = c(NA, NA, 0.578, 0.165, 0.246, 0.700, 0.348, 0.139, 1.017, NA, 0.531), 
#     lower = c(NA, NA, 0.372, 0.018, 0.072, 0.333, 0.083, 0.016, 0.365, NA, 0.386),
#     upper = c(NA, NA, 0.898, 1.517, 0.833, 1.474, 1.455, 1.209, 2.831, NA, 0.731)),
#     .Names = c("mean", "lower", "upper"), 
#     row.names = c(NA, -11L), 
#     class = "data.frame")
# 
# tabletext<-cbind(
#   c("", "Study", "Auckland", "Block", 
#     "Doran", "Gamsu", "Morrison", "Papageorgiou", 
#     "Tauesch", NA, "Summary"),
#   c("Deaths", "(steroid)", "36", "1", 
#     "4", "14", "3", "1", 
#     "8", NA, NA),
#   c("Deaths", "(placebo)", "60", "5", 
#     "11", "20", "7", "7", 
#     "10", NA, NA),
#   c("", "OR", "0.58", "0.16", 
#     "0.25", "0.70", "0.35", "0.14", 
#     "1.02", NA, "0.53"))
# 
# forestplot(tabletext, 
#            cochrane_from_rmeta,new_page = TRUE,
#            is.summary=c(TRUE,TRUE,rep(FALSE,8),TRUE),
#            clip=c(0.1,2.5), 
#            xlog=TRUE, 
#            col=fpColors(box="royalblue",line="darkblue", summary="royalblue"))