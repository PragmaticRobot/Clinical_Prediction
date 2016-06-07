### Let's see what PCA tells us! #########
## This script will pool the speed-related features and run PCA on them,
## see how much they covary


rm(list = ls()) # clear the environment

####################
## 99 - Functions ##
####################

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


#########################
## Setup and libraries ##
#########################

par.o <- par()
library(MASS)
library(R.matlab)
library(devtools)
library(rgl)
library(nFactors)
library(FactoMineR)
library(psych)
library(foreach)
library(randomForest)
library(doParallel)
library(inTrees)
library(lars)
library(glmnet)
library(coefplot)

##
setwd("C:/Users/Yaz/Documents/Clinical_Prediction/")
LASSO_res <- readMat("lasso_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
RF_res <- readMat("RF_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
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
# initialize data frames with the correct sizes and add row names for features
LASSO_1 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
LASSO_2 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
LASSO_3 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
RF_1 <- data.frame(matrix(vector(),55,100), row.names = gverb)
RF_2 <- data.frame(matrix(vector(),55,100), row.names = gverb)
RF_3 <- data.frame(matrix(vector(),55,100), row.names = gverb)


## put the features in the right order with coefficients (LASSO) 
# and variable importance (RF)
## after this, each matrix has features for rows, and cross-validation 
# results in columns
for (j in 1:100) # loop over cross-validations
{
  for (i in 1:dim(LASSO_res$LASSO.pred1[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred1[[2]][i,j])) {
      next 
    }
    else {
      LASSO_1[LASSO_res$LASSO.pred1[[2]][i,j],j] <- 
        LASSO_res$LASSO.pred1[[3]][i,j]
    }
  }
  ####################################################################
  for (i in 1:dim(LASSO_res$LASSO.pred2[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred2[[2]][i,j])) {
      next 
    }
    else {
      LASSO_2[LASSO_res$LASSO.pred2[[2]][i,j],j] <- 
        LASSO_res$LASSO.pred2[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:dim(LASSO_res$LASSO.pred3[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred3[[2]][i,j])) {
      next 
    }
    else {
      LASSO_3[LASSO_res$LASSO.pred3[[2]][i,j],j] <- 
        LASSO_res$LASSO.pred3[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:55)
  {
    RF_1[RF_res$RF.pred1[[2]][i,j],j] <- RF_res$RF.pred1[[3]][i,j]
    RF_2[RF_res$RF.pred2[[2]][i,j],j] <- RF_res$RF.pred2[[3]][i,j]
    RF_3[RF_res$RF.pred3[[2]][i,j],j] <- RF_res$RF.pred3[[3]][i,j]
  }
}

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

#########################################################
## The Section above creates two data frames that contain
## features for both the linear and quad models
## there is legacy code in that section, unimportant
## some things aren't being used, but no time to go through it
## Important outcomes from the above code:
## df1: linear data frame with subject x feature
## df2: quad data frame with subject x feature

## remove everything except df1 & 2 and variable names
# rm(list=ls()[! ls() %in% c("df1","df2","gverb","gm","yFM","yPartFM","yWO")])

## Current workspace:
## df1 & df2
## gm: variable names, abbreviated
## gverb: variable names, verbose

## speed related features in df1:
## 4,5,8,10 (spd pks), 17,18,22 (spd pks), 29,30,33,35
spd_cols <- df1[,c(4,5,8,10,17,18,22,29,30,33,35)]
spd.pca <- prcomp(spd_cols,center = TRUE, scale. = TRUE, rank. = 3)
summary(spd.pca) # variance accounted for
spd.pca$rotation # PC loadings
plot(spd.pca,type="lines") # scree plot
spd.pca$x # pca scores (values of observations for rotated data)
biplot(spd.pca) # plot pc1 against pc2
# ggbiplot(spd.pca, choices = 1:2,obs.scale = 1, var.scale = 1, 
         # ellipse = TRUE,circle = TRUE) # a prettier plot of pc1 vs pc2
plot3d(spd.pca$x[,1:3]) # a 3d plot of the pc1 - pc3

## factor analysis with varimax rotation
spd.fact <- factanal(spd_cols, 3, rotation = "varimax") # pick first 3 factors
print(spd.fact, digits = 2, cutoff = 0.3, sort = TRUE) # print the factors
load <- spd.fact$loadings[,1:2] # factor loadings of first two factors
plot(load, type = "p") # plot factor loadings
text(load, labels = names(spd_cols), cex = .7, pos = 3) # label features

## Determining the number of factors to extract
ev <- eigen(cor(spd_cols)) # get eigenvalues
ap <- parallel(subject=nrow(spd_cols),var=ncol(spd_cols),
               rep=100,cent=.05)
nS <- nScree(x=ev$values, aparallel=ap$eigen$qevpea)
par(par.o)
plotnScree(nS) 

## Using FactoMineR package
result_spd_pca <- PCA(spd_cols)

## whole data PCA
whole.pca <- prcomp(df1,center=TRUE, scale. = TRUE)
plot(whole.pca, type="l")
summary(whole.pca)
whole.pca$rotation
whole.pca$x # PCA scores
biplot(whole.pca)
plot3d(whole.pca$x[,1:3])
# How many factors?
ev0 <- eigen(cor(df1))
ap0 <- parallel(subject = nrow(df1), var = ncol(df1),
                rep=100,cent = 0.05)
ns0 <- nScree(x=ev0$values,aparallel = ap0$eigen$qevpea)
plotnScree(ns0)
# Varimax Rotation
### whole.fact <- factanal(df1, 4)
# since the matrix is singular, we need a different method, using the psych package
# options: fa, pa, minres, principal, factor.pa
# keep in mind we need 4 factors as per the previous step (plotnScree)
whole.fact <- fa(df1, nfactors = 4, residuals = TRUE, rotate = "varimax",
                            min.err = 0.001, fm="minres")
print(whole.fact$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)
op <- par(mar = c(4,20,4,2) + 0.1)
rownames(whole.fact$loadings) <- gverb
barplot(rev(whole.fact$loadings[,1]), horiz = TRUE, main="UEFM PC1", 
        names.arg =rev(rownames(whole.fact$loadings)),las=2)
barplot(rev(whole.fact$loadings[,2]), horiz = TRUE, main="UEFM PC2", 
        names.arg =rev(rownames(whole.fact$loadings)),las=2)
barplot(rev(whole.fact$loadings[,3]), horiz = TRUE, main="UEFM PC3", 
        names.arg =rev(rownames(whole.fact$loadings)),las=2)
barplot(rev(whole.fact$loadings[,4]), horiz = TRUE, main="UEFM PC4", 
        names.arg =rev(rownames(whole.fact$loadings)),las=2)

###########################################################
#### Now let's look at "rational" grouping of factors #####
###########################################################
par(par.o)
## Speed-related features: spd.fact
spd_cols <- df1[,c(4,5,8,10,17,18,22,29,30,33,35)]
stroke_cols <- df1[,c(42,44,45,46,47,48,49,50,51)] # stroke info
clinic_cols <- df1[,c(53,54,55)] # clinical Scores
demog_cols <- df1[,c(39,40,41,43)] # demographics
treat <- df1[,c(52)] # treatment
launch_cols <- df1[,c(1,3,7,14,16,20,26,28,32)]
mvmnt_cols <- df1[,c(2,6,9,13,15,19,21,25,27,31,34,38)]
accuracy_cols <- df1[,c(11,12,23,24,36,37)]

## The Seven Families:
# spd_cols
# stroke_cols
# clinic_cols
# demog_cols
# launch_cols
# mvmnt_cols
# accuracy_cols

## How many factors we need for each family?
ev1 <- eigen(cor(spd_cols)) # get eigenvalues
ev2 <- eigen(cor(stroke_cols)) # get eigenvalues
ev3 <- eigen(cor(clinic_cols)) # get eigenvalues
ev4 <- eigen(cor(demog_cols)) # get eigenvalues
ev5 <- eigen(cor(launch_cols)) # get eigenvalues
ev6 <- eigen(cor(mvmnt_cols)) # get eigenvalues
ev7 <- eigen(cor(accuracy_cols)) # get eigenvalues

ap1 <- parallel(subject=nrow(spd_cols),var=ncol(spd_cols), rep=100,cent=.05)
ap2 <- parallel(subject=nrow(stroke_cols),var=ncol(stroke_cols), rep=100,cent=.05)
ap3 <- parallel(subject=nrow(clinic_cols),var=ncol(clinic_cols), rep=100,cent=.05)
ap4 <- parallel(subject=nrow(demog_cols),var=ncol(demog_cols), rep=100,cent=.05)
ap5 <- parallel(subject=nrow(launch_cols),var=ncol(launch_cols), rep=100,cent=.05)
ap6 <- parallel(subject=nrow(mvmnt_cols),var=ncol(mvmnt_cols), rep=100,cent=.05)
ap7 <- parallel(subject=nrow(accuracy_cols),var=ncol(accuracy_cols), rep=100,cent=.05)

nS1 <- nScree(x=ev1$values, aparallel=ap1$eigen$qevpea)
nS2 <- nScree(x=ev2$values, aparallel=ap2$eigen$qevpea)
nS3 <- nScree(x=ev3$values, aparallel=ap3$eigen$qevpea)
nS4 <- nScree(x=ev4$values, aparallel=ap4$eigen$qevpea)
nS5 <- nScree(x=ev5$values, aparallel=ap5$eigen$qevpea)
nS6 <- nScree(x=ev6$values, aparallel=ap6$eigen$qevpea)
nS7 <- nScree(x=ev7$values, aparallel=ap7$eigen$qevpea)

plotnScree(nS1)
plotnScree(nS2)
plotnScree(nS3)
plotnScree(nS4)
plotnScree(nS5)
plotnScree(nS6)
plotnScree(nS7)

## Now we know how many factors:
# spd_cols 3
# stroke_cols 4
# clinic_cols 1
# demog_cols 1
# launch_cols 2
# mvmnt_cols 3
# accuracy_cols 2

## Let's do factor analysis with the requisite number of factors for each family:
spd_f <- fa(spd_cols, nfactors = 3, residuals = TRUE, rotate = "varimax", min.err = 0.001, fm="minres")
stroke_f <- fa(stroke_cols, nfactors = 4, residuals = TRUE, rotate = "varimax", min.err = 0.001, fm="minres")
clinic_f <- fa(clinic_cols, nfactors = 1, residuals = TRUE, rotate = "varimax", min.err = 0.001, fm="minres")
demog_f <- fa(demog_cols, nfactors = 1, residuals = TRUE, rotate = "varimax", min.err = 0.001, fm="minres")
launch_f <- fa(launch_cols, nfactors = 2, residuals = TRUE, rotate = "varimax", min.err = 0.001, fm="minres")
mvmnt_f <- fa(mvmnt_cols, nfactors = 3, residuals = TRUE, rotate = "varimax", min.err = 0.001, fm="minres")
accuracy_f <- fa(accuracy_cols, nfactors = 2, residuals = TRUE, rotate = "varimax", min.err = 0.001, fm="minres")

print(spd_f$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)
print(stroke_f$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)
print(clinic_f$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)
print(demog_f$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)
print(launch_f$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)
print(mvmnt_f$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)
print(accuracy_f$loadings, digits = 2, sort=FALSE, fixnames = TRUE, cutoff = 0.3)

spd_fa <- as.data.frame(spd_f$scores)
stroke_fa <- as.data.frame(stroke_f$scores)
clinic_fa <- as.data.frame(clinic_f$scores)
demog_fa <- as.data.frame(demog_f$scores)
launch_fa <- as.data.frame(launch_f$scores)
mvmnt_fa <- as.data.frame(mvmnt_f$scores)
accuracy_fa <- as.data.frame(accuracy_f$scores)
treat_fa <- as.data.frame(treat)

Factored_data <- cbind(spd_fa, stroke_fa, clinic_fa, demog_fa,
                       launch_fa, mvmnt_fa,accuracy_fa, treat_fa)
colnames(Factored_data) <- c("Speed Factor 3", 
                             "Speed Factor 1", 
                             "Speed Factor 2",
                             "Stroke Factor 4", 
                             "Stroke Factor 3", 
                             "Stroke Factor 2",
                             "Stroke Factor 1", 
                             "Clinical Scores Factor 1", 
                             "Demographics Factor 1",
                             "Launch Quality Factor 1", 
                             "Launch Quality Factor 2", 
                             "Whole Movement Factor 1",
                             "Whole Movement Factor 2", 
                             "Whole Movement Factor 3", 
                             "Accuracy Factor 1",
                             "Accuracy Factor 2",
                             "Treatment Flag")

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations

predLASSO1<-NULL
predLASSO2<-NULL
predLASSO3<-NULL

for (i in 1:100)
{
  predLASSO1[[i]] <- cv_mod(feat=Factored_data,out=yFM,stand=1,isWO = 0)
  predLASSO2[[i]] <- cv_mod(feat=Factored_data,out=yPartFM,stand=1,isWO = 0)
  predLASSO3[[i]] <- cv_mod(feat=Factored_data,out=yWO,stand=1,isWO = 1)
}
LASSO_pred1 <- cleanup_LASSO(predLASSO1)
LASSO_pred2 <- cleanup_LASSO(predLASSO2)
LASSO_pred3 <- cleanup_LASSO(predLASSO3)

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
