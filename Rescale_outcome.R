# This script will rescale the outcome (change in clinical score) to try to account for
# non-normality as per Konrad's suggestion. This reconstructs the models, then plots the new
# fits and histograms.

# Two transformations of the outcome are suggested: relative change in score
# and log(change in score)
# just because, will not do Arm-only FM

#########################
## Setup and libraries ##
#########################

rm(list = ls())

####################################
######### 99 - Functions ###########
####################################

cv_RF <- function(feat,out)
{
  RF_model <- foreach(ntree = rep(10000,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% randomForest(feat,out,ntree=ntree)
  imp <- sort(RF_model$importance, partial=NULL,decreasing=TRUE, index.return=TRUE)
  predict <- predict(RF_model, feat, type = "response" , predict.all=FALSE)
  GRE <- list(predict=predict, Index=imp$ix, imps=imp$x)
  return(GRE)
}

cv_mod <- function(feat, out, stand, isWO)
{
  if(stand == 1)
  {
    lassoMOD <- cv.glmnet(x=as.matrix(feat),y=out, standardize = TRUE)
  } else {
    lassoMOD <- cv.glmnet(x=as.matrix(feat),y=out, standardize = FALSE)
  }
  cutoff <- min(lassoMOD$cvup)
  means <- lassoMOD$cvm
  if(isWO == 1)
  {
    ind <- min(which(means<cutoff)) # now we know which model won
  } else {
    ind <- max(which(means<cutoff)) # now we know which model won
  }
  predict <- predict(lassoMOD, newx=as.matrix(feat), s= lassoMOD$lambda[ind])
  cc <- coef(lassoMOD$glmnet.fit, s = lassoMOD$lambda[ind])
  summ <- summary(cc)
  GRE <- list(predict=predict, Index = summ$i, coefs = summ$x)
  return(GRE)
}

## for the cleanup functions flag:
## 1: predictions
## 2: indices of winning variables (LASSO) or order of variable importance (RF)
## 3: value for coefficients of winning variables, or importance of variable (RF)
cleanup_LASSO <- function(pred1)
{
  predictions1 <- pred1[[1]]$predict
  inds1 <- pred1[[1]]$Index
  coefs1 <- pred1[[1]]$coefs
  for (i in 2:100)
  {
    predictions1 <- cbind(predictions1,pred1[[i]]$predict)
    inds1 <- qpcR:::cbind.na(inds1,pred1[[i]]$Index)
    coefs1 <- qpcR:::cbind.na(coefs1,pred1[[i]]$coefs)
  }
  LASSO_pred1 <- list(Predictions=predictions1,Indices=inds1,Coeffs=coefs1)
  return(LASSO_pred1)
}

cleanup_RF <- function(pred1)
{
  predictions1 <- pred1[[1]]$predict
  inds1 <- pred1[[1]]$Index
  imps1 <- pred1[[1]]$imps
  for (i in 2:100)
  {
    predictions1 <- cbind(predictions1,pred1[[i]]$predict)
    inds1 <- qpcR:::cbind.na(inds1,pred1[[i]]$Index)
    imps1 <- qpcR:::cbind.na(imps1,pred1[[i]]$imps)
  }
  RF_pred1 <- list(Predictions=predictions1,Indices=inds1,Imps=imps1)
  return(RF_pred1)
}

# plot_sideways plots the variable importance for the cross-validation of
# the random forest algorithm, it takes the "sorted list", a data frame with
# representing features and columns representing cross-validations, the rows
# should have the feature name as rownames
#
# The second input is N, if N is zero, the function uses the number of features
# in the list as N, otherwise it uses the value provided when the function is
# called
#
# The third input is a title string, used for the "main" title for the plot
#
#
plot_sideways <- function(sorted_list,N,tit)
{
  par(las=1)
  par(mar=c(4,20,4,2))
  x <- sorted_list[1,]
  nam <- rownames(sorted_list)[1]
  x <- x[!is.na(x)]
  offs <- runif(length(x),-0.35,0.35)
  gg <- x
  if(N == 0){
    N <- dim(sorted_list)[1]
  }
  x1 <- ceiling(max(sorted_list, na.rm = TRUE))
  x2 <- floor(min(sorted_list, na.rm = TRUE))
  yy <- rep(N, length(gg))+ offs
  plot(gg,yy,col='blue',pch=16,ylim = c(1, N),
       xlim = c(x2, x1),cex=0.3,yaxt='n', main=tit,
       xlab="",ylab="",cex.lab=0.7)
  # mean segment
  xx1 <- fivenum(gg,na.rm = TRUE)
  xq1 <- xx1[2];  xq2 <- xx1[3]; xq3 <- xx1[4]
  segments(xq2-0.05, N-0.5, xq2-0.05, N+0.5, col = 'red',lwd = 2) # center (median)
  segments(xq1-0.05, N-0.5, xq1-0.05, N+0.5, col = 'forestgreen',lwd = 2) # left (Q1)
  segments(xq3-0.05, N-0.5, xq3-0.05, N+0.5, col = 'forestgreen',lwd = 2) # right (Q3)
  segments(xq1-0.05, N+0.5, xq3-0.05, N+0.5, col = 'forestgreen',lwd = 2) # top
  segments(xq1-0.05, N-0.5, xq3-0.05, N-0.5, col = 'forestgreen',lwd = 2) # bottom
  axis(side = 2, at = N,paste(nam))
  for (i in 2:N){
    x <- sorted_list[i,]
    nam <- rownames(sorted_list)[i]
    x <- x[!is.na(x)]
    offs <- runif(length(x),-0.35,0.35)
    gg <- x
    yy <- rep(N - (i-1), length(gg))+ offs
    points(gg,yy,col='blue',pch=16,cex=0.3)
    yc <- N - (i-1) # the center of the y-axis for this feature
    xx1 <- fivenum(gg,na.rm = TRUE)
    xq1 <- xx1[2];  xq2 <- xx1[3]; xq3 <- xx1[4]
    segments(xq2-0.05, yc-0.5, xq2-0.05, yc+0.5, col = 'red',lwd = 2) # center (median)
    segments(xq1-0.05, yc-0.5, xq1-0.05, yc+0.5, col = 'forestgreen',lwd = 2) # left (Q1)
    segments(xq3-0.05, yc-0.5, xq3-0.05, yc+0.5, col = 'forestgreen',lwd = 2) # right (Q3)
    segments(xq1-0.05, yc+0.5, xq3-0.05, yc+0.5, col = 'forestgreen',lwd = 2) # top
    segments(xq1-0.05, yc-0.5, xq3-0.05, yc-0.5, col = 'forestgreen',lwd = 2) # bottom
    axis(side = 2, at = yc,paste(nam))
  }
}

#################################
########## Libraries ############
#################################
# options(error=recover)
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,
               psych, foreach, randomForest, doParallel, inTrees,tableone,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2)

###############################################
############ Create the data frames ###########
###############################################
##
load("colnames.rda")
## restructure the data into a data frame
# RF: Random Forests
# LS: LASSO
################################## IMPORTANT NOTE ########
#### This section ONLY deals with LINEAR results, 
#### quadratic models are not included!

gm2 <- c("Intercept",gm)    # add "intercept" to the list of feature names for LASSO
gverb <- gm                 # verbose version of feature names

######### Making feature names verbose (looong list) #############
gverb[1] <- "Mean Reaction Time"; gverb[2] <- "Mean Trial Time"
gverb[3] <- "Mean Initial Direction Error"; gverb[4] <- "Mean Pre-Movement Speed"
gverb[5] <- "Mean Max Speed"; gverb[6] <- "Mean Hand Path Length"
gverb[7] <- "Mean Initial Movement Ratio"; gverb[8] <- "Mean Speed Ratio"
gverb[9] <- "Mean Path Length Ratio"; gverb[10] <- "Mean Number of Speed Peaks"
gverb[11] <- "Mean Max Perpendicular Distance"; 
gverb[12] <- "Mean Percent of Movement in the Target Direction"
gverb[13] <- "Mean Mean Arrest Period Ratio"; gverb[14] <- "Max Reaction Time"
gverb[15] <- "Max Trial Time"; gverb[16] <- "Max Initial Direction Error"
gverb[17] <- "Max Pre-Movement Speed"; gverb[18] <- "Max Max Speed"
gverb[19] <- "Max Hand Path Length"; gverb[20] <- "Max Initial Movement Ratio"
gverb[21] <- "Max Path Length Ratio"; gverb[22] <- "Max Number of Speed Peaks"
gverb[23] <- "Max Max Perpendicular Distance"; 
gverb[24] <- "Max Percent of Movement in the Target Direction"
gverb[25] <- "Max Mean Arrest Period Ratio"; gverb[26] <- "Var Reaction Time"
gverb[27] <- "Var Trial Time"; gverb[28] <- "Var Initial Direction Error"
gverb[29] <- "Var Pre-Movement Speed"; gverb[30] <- "Var Max Speed"
gverb[31] <- "Var Hand Path Length"; gverb[32] <- "Var Initial Movement Ratio"
gverb[33] <- "Var Speed Ratio"; gverb[34] <- "Var Path Length Ratio"
gverb[35] <- "Var Number of Speed Peaks"; 
gverb[36] <- "Var Max Perpendicular Distance"
gverb[37] <- "Var Percent of Movement in the Target Direction"; 
gverb[38] <- "Var Mean Arrest Period Ratio"
gverb[39] <- "Age"; gverb[40] <- "Height"; gverb[41] <- "Mass"; 
gverb[42] <- "Months Post-Stroke"; gverb[43] <- "Sex"; 
gverb[44] <- "Left Hand Dominant?"
gverb[45] <- "Left Side Affected?"; gverb[46] <- "Dominant Side Affected?"
gverb[47] <- "Hemorrhagic Stroke?"; gverb[48] <- "Unknown Stroke Location?"
gverb[49] <- "Cortical Stroke?"; gverb[50] <- "Subcortical Stroke?"
gverb[51] <- "Brainstem Stroke?"; gverb[52] <- "Error-Augmentation Treatment?"
gverb[53] <- "Initial Fugl-Meyer Score"; 
gverb[54] <- "Initial Wolf Motor Function Score"
gverb[55] <- "Initial Box-and-Blocks Score"

##
gverb2 <- c("Intercept",gverb) # add intercept

#####################################################
## New Wroking folder: Dropbox/Research/NewDec2015 ##
##                 Data Prep - Importing           ##
#####################################################

Data <- readMat("RF_Ready_noscramble.mat", maxLength=NULL, fixNames=TRUE, 
                verbose=FALSE)
# setup variables
Features <- data.frame(Data$red.DesignMat[1:26,2:56])
Outcomes <- data.frame(Data$y.FM,Data$y.PartsFM,Data$y.WO)
numfeat <- dim(Features)[2]

###############################
## Data Prep - column titles ##
###############################

colnames(Features)[1] <- Data$red.Dnames.full[[2]][[1]][[1,1]]
for (i in 2:numfeat+1)
{
  colnames(Features)[i-1] <- Data$red.Dnames.full[[i]][[1]][[1,1]]
}
colnames(Outcomes)[1] <- 'FullFM'
colnames(Outcomes)[2] <- 'ArmFM'
CombinedNum <- cbind(Features[,1:42],Features[,53:55],Outcomes)
CombinedFull <- cbind(Features,Outcomes)
rm(Data,Dpath,i,numfeat)

colnames(Features) <- gsub(" ","",colnames(Features)) # ggpairs 
# doesn't like spaces in varnames
top <- dim(Features)[2]

## We want to use the "Features" Data Frame for the next part(s)

############################################################
## Invert binary variables (there was an error in matlab) ##
############################################################
VarNames <- colnames(Features)
varsToFactor <- VarNames[43:52]
FeaturesF <- Features # a copy of features with binary variables 
# turned into factors (so the original is fine)
FeaturesF[varsToFactor] <- lapply(FeaturesF[varsToFactor], factor)

for (i in 43:52)
{
  FeaturesF[,i] <- as.factor(2-as.integer(FeaturesF[,i])) 
}

outs <- list()
outs[[1]] <- Outcomes$FullFM
yFM = Outcomes$FullFM
outs[[2]] <- Outcomes$ArmFM
yPartFM = Outcomes$ArmFM
outs[[3]] <- Outcomes$Data.y.WO
yWO = Outcomes$Data.y.WO

###############################################
## 03 - Prepare for quadratic random forests ##
###############################################
df1 <- Features 
# this is the original copy of Features, where binary variables 
# are still treated as integers
# this allows for the multiplication of binary variables with other variables

df2 <- Features
names2 <- c(colnames(Features))
for (i in c(1:dim(df1)[2])) {
  for (j in c(i:dim(df1)[2])){
    name_new <- paste(colnames(df1)[i], '*',colnames(df1)[j], sep="")
    names2 <- c(names2,name_new)
    
    df_new <- data.frame(new_name = t(t(df1[,i])*t(df1[,j])))
    df2 <- cbind(df2,(df_new))
  }
}
colnames(df2) <- c(names2)
FeaturesQ <- df2               
# a new data frame containing quadratic predictors with correct column
# names (thanks Peter!)
FeaturesQ[varsToFactor] <- lapply(FeaturesQ[varsToFactor], factor)

for (i in 43:52)
{
  FeaturesQ[,i] <- as.factor(2-as.integer(FeaturesQ[,i])) 
}

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations

lgFM <- log(yFM)
lgWO <- log(yWO)

rvFM <- yFM/df1$init_FM
rvWO <- yWO/df1$initWO

predLASSO_lg1<-NULL
predLASSO_lg2<-NULL
predLASSO_lg3<-NULL
predLASSO_lg4<-NULL

predLASSO_rv1<-NULL
predLASSO_rv2<-NULL
predLASSO_rv3<-NULL
predLASSO_rv4<-NULL

for (i in 1:100)
{
  predLASSO_lg1[[i]] <- cv_mod(feat=df1,out=lgFM,stand=1,isWO = 0)
  predLASSO_lg2[[i]] <- cv_mod(feat=df1,out=lgWO,stand=1,isWO = 1)
  predLASSO_lg3[[i]] <- cv_mod(feat=df2,out=lgFM,stand=1,isWO = 0)
  predLASSO_lg4[[i]] <- cv_mod(feat=df2,out=lgWO,stand=1,isWO = 1)
  predLASSO_rv1[[i]] <- cv_mod(feat=df1,out=rvFM,stand=1,isWO = 0)
  predLASSO_rv2[[i]] <- cv_mod(feat=df1,out=rvWO,stand=1,isWO = 1)
  predLASSO_rv3[[i]] <- cv_mod(feat=df2,out=rvFM,stand=1,isWO = 0)
  predLASSO_rv4[[i]] <- cv_mod(feat=df2,out=rvWO,stand=1,isWO = 1)
}

LASSO_pred_lg1 <- cleanup_LASSO(predLASSO_lg1)
LASSO_pred_lg2 <- cleanup_LASSO(predLASSO_lg2)
LASSO_pred_lg3 <- cleanup_LASSO(predLASSO_lg3)
LASSO_pred_lg4 <- cleanup_LASSO(predLASSO_lg4)
LASSO_pred_rv1 <- cleanup_LASSO(predLASSO_rv1)
LASSO_pred_rv2 <- cleanup_LASSO(predLASSO_rv2)
LASSO_pred_rv3 <- cleanup_LASSO(predLASSO_rv3)
LASSO_pred_rv4 <- cleanup_LASSO(predLASSO_rv4)

predRF1 <- NULL
predRF2 <- NULL
predRF3 <- NULL
for (i in 1:100)
{
  predRF1[[i]] <- cv_RF(feat=Factored_data, out=yFM)
  predRF2[[i]] <- cv_RF(feat=Factored_data, out=yPartFM)
  predRF3[[i]] <- cv_RF(feat=Factored_data, out=yWO)
}

RF_pred1 <- cleanup_RF(predRF1)
RF_pred2 <- cleanup_RF(predRF2)
RF_pred3 <- cleanup_RF(predRF3)

## LASSO
FM_lin_lasso <-   LASSO_pred1[[1]]
pFM_lin_lasso <-  LASSO_pred2[[1]]
WO_lin_lasso  <-  LASSO_pred3[[1]]

## Random Forests
FM_lin_RF <-   RF_pred1[[1]]
pFM_lin_RF <-  RF_pred2[[1]]
WO_lin_RF <-   RF_pred3[[1]]
