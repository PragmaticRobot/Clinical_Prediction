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
library(qpcR)
library(dplyr)
library(lattice) 

##
Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","lasso_pred.mat")
LASSO_res <- readMat(Dpath, maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
Dpath <- file.path("C:","Users","Yaz","Dropbox","Research","NewDec2015","2016-04-19 - Second Run CV","RF_pred.mat")
RF_res <- readMat(Dpath, maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
load(file = "C:/Users/Yaz/Dropbox/Research/NewDec2015/colnames.rda")
## restructure the data into a data frame
# RF: Random Forests
# LS: LASSO
################################## IMPORTANT NOTE #################################################
######### This section ONLY deals with LINEAR results, quadratic models are not included! #########
###################################################################################################
gm2 <- c("Intercept",gm)                                        # add "intercept" to the list of feature names
# initialize data frames with the correct sizes and add row names for features
LASSO_1 <- data.frame(matrix(vector(),56,100), row.names = gm2)
LASSO_2 <- data.frame(matrix(vector(),56,100), row.names = gm2)
LASSO_3 <- data.frame(matrix(vector(),56,100), row.names = gm2)
RF_1 <- data.frame(matrix(vector(),55,100), row.names = gm)
RF_2 <- data.frame(matrix(vector(),55,100), row.names = gm)
RF_3 <- data.frame(matrix(vector(),55,100), row.names = gm)


## put the features in the right order with coefficients (LASSO) and variable importance (RF)
## after this, each matrix has features for rows, and cross-validation results in columns
for (j in 1:100) # loop over cross-validations
{
  for (i in 1:dim(LASSO_res$LASSO.pred1[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred1[[2]][i,j])) {
      next 
    }
    else {
      LASSO_1[LASSO_res$LASSO.pred1[[2]][i,j],j] <- LASSO_res$LASSO.pred1[[3]][i,j]
    }
  }
  ####################################################################
  for (i in 1:dim(LASSO_res$LASSO.pred2[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred2[[2]][i,j])) {
      next 
    }
    else {
      LASSO_2[LASSO_res$LASSO.pred2[[2]][i,j],j] <- LASSO_res$LASSO.pred2[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:dim(LASSO_res$LASSO.pred3[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred3[[2]][i,j])) {
      next 
    }
    else {
      LASSO_3[LASSO_res$LASSO.pred3[[2]][i,j],j] <- LASSO_res$LASSO.pred3[[3]][i,j]
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
counts_FM <- data.frame(counts_FM[2:length(counts_FM)], row.names = gm)
colnames(counts_FM) <- c("Count")
counts_pFM <- data.frame(counts_pFM[2:length(counts_pFM)], row.names = gm)
colnames(counts_pFM) <- c("Count")
counts_WO <- data.frame(counts_WO[2:length(counts_WO)], row.names = gm)
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
op <- par(mar = c(4,10,4,2) + 0.1)

# plot bar chart, rev makes sure the first item is on top (reverse, since default is botton)
barplot(rev(sorted_FM$Count), horiz = TRUE, main="FM Variable Frequency", names.arg =rev(rownames(sorted_FM)),las=2)
barplot(rev(sorted_pFM$Count), horiz = TRUE, main="Arm-Only FM Variable Frequency", names.arg =rev(rownames(sorted_pFM)),las=2)
barplot(rev(sorted_WO$Count), horiz = TRUE, main="Wolf Variable Frequency", names.arg =rev(rownames(sorted_WO)),las=2)
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
FM_Lin_LASSO_sort <- gg[order(-abs(gg$X1)), , drop = FALSE]

## now repeat for Arm FM
g<-which(counts_pFM==0,arr.ind = T)
gnames <- namelist[-c(as.vector(g[,1]))]
LASSO_2 <- data.frame(LASSO_2[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_2[-c(as.vector(g[,1])),],row.names = gnames)
ArmFM_Lin_LASSO_sort <- gg[order(-abs(gg$X1)), , drop = FALSE]

## and Wolf
g<-which(counts_WO==0,arr.ind = T)
gnames <- namelist[-c(as.vector(g[,1]))]
LASSO_3 <- data.frame(LASSO_3[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_3[-c(as.vector(g[,1])),],row.names = gnames)
Wolf_Lin_LASSO_sort <- gg[order(-abs(gg$X1)), , drop = FALSE]

## NOW, sort RF results, easier since no empty cells
FM_Lin_RF_sort <- RF_1[order(-RF_1$X1), , drop=FALSE]
ArmFM_Lin_RF_sort <- RF_2[order(-RF_2$X1), , drop=FALSE]
Wolf_Lin_RF_sort <- RF_3[order(-RF_3$X1), , drop=FALSE]
