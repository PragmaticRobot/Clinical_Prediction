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
library(qpcR)

##
setwd("C:/Users/Yaz/Documents/Clinical_Prediction/")
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

## LASSO
FM_lin_lasso <-   LASSO_pred1[[1]]
pFM_lin_lasso <-  LASSO_pred2[[1]]
WO_lin_lasso  <-  LASSO_pred3[[1]]

## Random Forests
FM_lin_RF <-   RF_pred1[[1]]
pFM_lin_RF <-  RF_pred2[[1]]
WO_lin_RF <-   RF_pred3[[1]]

## fitting the linear models
fm_lin_lasso_mod <- lm(yFM~rowMeans(FM_lin_lasso))
fm_lin_rf_mod <- lm(yFM~rowMeans(FM_lin_RF))

pfm_lin_lasso_mod <- lm(yPartFM~rowMeans(pFM_lin_lasso))
pfm_lin_rf_mod <- lm(yPartFM~rowMeans(pFM_lin_RF))

wo_lin_lasso_mod <- lm(yWO~rowMeans(WO_lin_lasso))
wo_lin_rf_mod <- lm(yWO~rowMeans(WO_lin_RF))

## plotting I: Fugl-Meyer, make plot pretty
plot(rowMeans(FM_lin_lasso),yFM, ylim = c(-10, 20), xlim = c(-10, 20), col='red',
     xlab="Predicted Change in Fugl-Meyer Score",ylab="Change in Fugl-Meyer Score",
     cex.lab=0.8,pch=21, cex=0.8)
abline(fm_lin_lasso_mod, col='red')
blah <- summary(fm_lin_lasso_mod)
r1 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(FM_lin_RF),yFM, col = 'blue',pch=24, cex=0.8)
abline(fm_lin_rf_mod, col='blue')
blah <- summary(fm_lin_rf_mod)
r3 <- round(blah$adj.r.squared, digits = 2)

# abline(0,1, col='black',pch=16,cex=0.6)
legend("topleft",c(paste('Linear LASSO R^2 = ',as.character(r1)),
                   paste('Linear Random Forests R^2 = ',as.character(r3))),
       cex=0.8, col=c('red','blue'), text.col = 'black',
       pch=c(21,24),bty="n")


## plotting II: Arm-only Fugl-Meyer, make plot pretty
plot(rowMeans(pFM_lin_lasso),yPartFM, ylim = c(-7, 12), xlim = c(-7, 12), col='red',
     xlab="Predicted Change in Arm Fugl-Meyer Score",ylab="Change in Arm-Only Fugl-Meyer Score",
     cex.lab=0.8,pch=21, cex=0.8)
abline(pfm_lin_lasso_mod, col='red')
blah <- summary(pfm_lin_lasso_mod)
r5 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(pFM_lin_RF),yPartFM, col = 'blue',pch=24, cex=0.8)
abline(pfm_lin_rf_mod, col='blue')
blah <- summary(pfm_lin_rf_mod)
r7 <- round(blah$adj.r.squared, digits = 2)

# abline(0,1, col='black',pch=16,cex=0.6)
legend("topleft",c(paste('Linear LASSO R^2 = ',as.character(r5)),
                   paste('Linear Random Forests R^2 = ',as.character(r7))),
       cex=0.8, col=c('red','blue'), text.col = 'black',
       pch=c(21,24),bty="n")


## plotting III: Wolf, make plot pretty
plot(rowMeans(WO_lin_lasso),yWO, ylim = c(-28, 14), xlim = c(-28, 14), col='red',
     xlab="Predicted Change in Wolf Score",ylab="Change in Wolf Motor Function Score",
     cex.lab=0.8,pch=21, cex=0.8)
abline(wo_lin_lasso_mod, col='red')
blah <- summary(wo_lin_lasso_mod)
r9 <- round(blah$adj.r.squared, digits = 2)
points(rowMeans(WO_lin_RF),yWO, col = 'blue',pch=24, cex=0.8)
abline(wo_lin_rf_mod, col='blue')
blah <- summary(wo_lin_rf_mod)
r11 <- round(blah$adj.r.squared, digits = 2)

# abline(0,1, col='black',pch=16,cex=0.6)
legend("topleft",c(paste('Linear LASSO R^2 = ',as.character(r9)),
                   paste('Linear Random Forests R^2 = ',as.character(r11))),
       cex=0.8, col=c('red','blue'), text.col = 'black',
       pch=c(21,24),bty="n")


coln1 <- colnames(Factored_data)
coln2 <- c("Intercept",coln1)
# initialize data frames with the correct sizes and add row names for features
LASSO_1 <- data.frame(matrix(vector(),18,100), row.names = coln2)
LASSO_2 <- data.frame(matrix(vector(),18,100), row.names = coln2)
LASSO_3 <- data.frame(matrix(vector(),18,100), row.names = coln2)
RF_1 <- data.frame(matrix(vector(),17,100), row.names = coln1)
RF_2 <- data.frame(matrix(vector(),17,100), row.names = coln1)
RF_3 <- data.frame(matrix(vector(),17,100), row.names = coln1)


## put the features in the right order with coefficients (LASSO) 
# and variable importance (RF)
## after this, each matrix has features for rows, and cross-validation 
# results in columns
for (j in 1:100) # loop over cross-validations
{
  for (i in 1:dim(LASSO_pred1[[2]])[1])
  {
    if (is.na(LASSO_pred1[[2]][i,j])) {
      next 
    }
    else {
      LASSO_1[LASSO_pred1[[2]][i,j],j] <- 
        LASSO_pred1[[3]][i,j]
    }
  }
  ####################################################################
  for (i in 1:dim(LASSO_pred2[[2]])[1])
  {
    if (is.na(LASSO_pred2[[2]][i,j])) {
      next 
    }
    else {
      LASSO_2[LASSO_pred2[[2]][i,j],j] <- 
        LASSO_pred2[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:dim(LASSO_pred3[[2]])[1])
  {
    if (is.na(LASSO_pred3[[2]][i,j])) {
      next 
    }
    else {
      LASSO_3[LASSO_pred3[[2]][i,j],j] <- 
        LASSO_pred3[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:17)
  {
    RF_1[RF_pred1[[2]][i,j],j] <- RF_pred1[[3]][i,j]
    RF_2[RF_pred2[[2]][i,j],j] <- RF_pred2[[3]][i,j]
    RF_3[RF_pred3[[2]][i,j],j] <- RF_pred3[[3]][i,j]
  }
}

## Calculating rank based on frequency of variable selection for LASSO

# Initialize the count holders
counts_FM <- rep(0, dim(LASSO_3)[[1]])
counts_pFM <- rep(0, dim(LASSO_3)[[1]])
counts_WO <- rep(0, dim(LASSO_3)[[1]])

# count
for (ii in 1:length(counts_FM)){
  for (j in 1:100){
    if (is.na(LASSO_1[[ii,j]]) || LASSO_1[[ii,j]]==0){
      counts_FM[[ii]] <- counts_FM[[ii]] + 0
    } else {
      counts_FM[[ii]] <- counts_FM[[ii]]+1
    }
    if (is.na(LASSO_2[[ii,j]]) || LASSO_2[[ii,j]]==0){
      counts_pFM[[ii]] <- counts_pFM[[ii]] + 0
    } else {
      counts_pFM[[ii]] <- counts_pFM[[ii]]+1
    }
    if (is.na(LASSO_3[[ii,j]]) || LASSO_3[[ii,j]]==0){
      counts_WO[[ii]] <- counts_WO[[ii]] + 0
    } else {
      counts_WO[[ii]] <- counts_WO[[ii]]+1
    }
  }
}

# convert the count holders to data frams (to have row names)
counts_FM <- data.frame(counts_FM[2:length(counts_FM)], row.names = coln1)
colnames(counts_FM) <- c("Count")
counts_pFM <- data.frame(counts_pFM[2:length(counts_pFM)], row.names = coln1)
colnames(counts_pFM) <- c("Count")
counts_WO <- data.frame(counts_WO[2:length(counts_WO)], row.names = coln1)
colnames(counts_WO) <- c("Count")

# sort and pick out non-zero elements
sorted_FM<- counts_FM[order(-counts_FM$Count), , drop = FALSE]
sorted_pFM<- counts_pFM[order(-counts_pFM$Count), , drop = FALSE]
sorted_WO<- counts_WO[order(-counts_WO$Count), , drop = FALSE]
sorted_FM<- sorted_FM[seq(17), , FALSE]
sorted_pFM<-sorted_pFM[seq(9), , FALSE]
sorted_WO<- sorted_WO[seq(23), , FALSE]

# Simple Bar Plot with Added Labels

# set the plot margins to have more space on the left for variable names
op <- par(mar = c(4,20,4,2) + 0.1)

# plot bar chart, rev makes sure the first item is on top (reverse, since default is botton)
# png("UEFM_counts.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
barplot(rev(sorted_FM$Count), horiz = TRUE, main="UEFM Variable Frequency", names.arg =rev(rownames(sorted_FM)),las=2)
# dev.off()

# png("ArmFM_counts.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
barplot(rev(sorted_pFM$Count), horiz = TRUE, main="Arm-Only FM Variable Frequency", names.arg =rev(rownames(sorted_pFM)),las=2)
# dev.off()

# png("Wolf_counts.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
barplot(rev(sorted_WO$Count), horiz = TRUE, main="Wolf Variable Frequency", names.arg =rev(rownames(sorted_WO)),las=2)
# dev.off()

par(op) ## reset plot parameters (margins)

## sort LASSO results (use first column for sorting)
# first, find the rows that have at leasst one value
namelist <- rownames(counts_FM)                       # full list of variable names
# g will hold indices of zero frequency variables
g<-which(counts_FM==0,arr.ind = T)
# gnames will contain the names of the "remaining" features (non-zero freq)
gnames <- namelist[-c(as.vector(g[,1]))]
# remove the first row (intercept) from all three LASSO matrices
LASSO_1 <- data.frame(LASSO_1[-c(1),],row.names = namelist)
# gg will be the data frame without the unwanted features
gg <- data.frame(LASSO_1[-c(as.vector(g[,1])),],row.names = gnames)
mns <- rowMeans(gg, na.rm=TRUE);
FM_Lin_LASSO_sort <- gg[order(-abs(mns)), , drop = FALSE]

## now repeat for Arm FM
g<-which(counts_pFM==0,arr.ind = T)
gnames <- namelist[-c(as.vector(g[,1]))]
LASSO_2 <- data.frame(LASSO_2[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_2[-c(as.vector(g[,1])),],row.names = gnames)
mns <- rowMeans(gg, na.rm=TRUE);
ArmFM_Lin_LASSO_sort <- gg[order(-abs(mns)), , drop = FALSE]

## and Wolf
g<-which(counts_WO==0,arr.ind = T)
gnames <- namelist[-c(as.vector(g[,1]))]
LASSO_3 <- data.frame(LASSO_3[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_3[-c(as.vector(g[,1])),],row.names = gnames)
mns <- rowMeans(gg, na.rm=TRUE);
Wolf_Lin_LASSO_sort <- gg[order(-abs(mns)), , drop = FALSE]

## NOW, sort RF results, easier since no empty cells
FM_Lin_RF_sort <- RF_1[order(-RF_1$X1), , drop=FALSE]
ArmFM_Lin_RF_sort <- RF_2[order(-RF_2$X1), , drop=FALSE]
Wolf_Lin_RF_sort <- RF_3[order(-RF_3$X1), , drop=FALSE]

#################### plot the results #################################

## Fugl-Meyer:
tit1 <- "UEFM LASSO Coefficients"
tit2 <- "UEFM RF Variable Importance"

# png("UEFM_LASSO.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(FM_Lin_LASSO_sort,0,tit1)
# dev.off()

# png("UEFM_RF.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(FM_Lin_RF_sort,17,tit2)
# dev.off()

## Arm-Only FM:
tit1 <- "Arm-Only FM LASSO Coefficients"
tit2 <- "Arm-Only FM RF Variable Importance"
# png("ArmFM_LASSO.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(ArmFM_Lin_LASSO_sort,0,tit1)
# dev.off()

# png("ArmFM_RF.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(ArmFM_Lin_RF_sort,17,tit2)
# dev.off()

## Wolf:
tit1 <- "Wolf LASSO Coefficients"
tit2 <- "Wolf RF Variable Importance"
# png("Wolf_LASSO.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(Wolf_Lin_LASSO_sort,0,tit1)
# dev.off()

# png("Wolf_RF.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(Wolf_Lin_RF_sort,17,tit2)
# dev.off()
