# This version adds bootstrapping to LASSO, so the model predictions
# can be averaged and mean and std can be calculated for these
# predictions, the results will be forwarded to matlab for plotting

rm(list = ls()) # clear the environment

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
# source("functions2.R")
load("2016dec12_Basics.rda")
source('create_df1_df2_no_unknown.R')

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations

####################################
# MAT Files:
# LASSO_res <- readMat("lasso_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
# RF_res <- readMat("RF_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
# RDA Files
load("2016dec12_RF_nReps100_nFolds4.rda")
load("2016dec12_LS_nReps100_nFolds4.rda")
# load("colnames.rda")
## restructure the data into a data frame
# RF: Random Forests
# LS: LASSO
# initialize data frames with the correct sizes and add row names for features
colnmquad <- c("Intercept",colnames(df2)) #colnames for quadratic models
colnmNoExtrap <- c("Intercept",colnames(df3))
# gverb2 <- gverb2[-c(49)] # remove unknown stroke location feature
# gverb <- gverb[-c(48)]
LASSO_1 <- data.frame(matrix(vector(),52,100), row.names = gverb2)
LASSO_2 <- data.frame(matrix(vector(),52,100), row.names = gverb2)
LASSO_5 <- data.frame(matrix(vector(),1366,100),row.names = colnmNoExtrap)
LASSO_6 <- data.frame(matrix(vector(),1366,100),row.names = colnmNoExtrap)

RF_3 <- data.frame(matrix(vector(),51,100), row.names = gverb)
RF_4 <- data.frame(matrix(vector(),51,100), row.names = gverb)
RF_7 <- data.frame(matrix(vector(),1365,100), row.names = colnames(df3))
RF_8 <- data.frame(matrix(vector(),1365,100), row.names = colnames(df3))

## put the features in the right order with coefficients (LASSO) and variable importance (RF)
## after this, each matrix has features for rows, and cross-validation results in columns
# LASSO_1 <- reorder_LS(LASSO_pred1)
# LASSO_2 <- reorder_LS(LASSO_pred5)
# LASSO_5 <- reorder_LS(LASSO_pred2)
# LASSO_6 <- reorder_LS(LASSO_pred6)
# 
# RF_3 <- reorder_RF(RF_pred3)
# RF_4 <- reorder_RF(RF_pred7)
# RF_7 <- reorder_RF(RF_pred4)
# RF_8 <- reorder_RF(RF_pred8)

#############################
##### Fix the reorder #######
#############################
## 2 and 5 are reversed because we want LS_L WO_L LS_Q WO_Q 
## (different order from Step1_models.R)

for (j in 1:100) # loop over cross-validations
{
  # for (i in 1:dim(LASSO_res$LASSO.pred1[[2]])[1])
  for (i in 1:dim(LASSO_pred1[[2]])[1])
  {
    # if (is.na(LASSO_res$LASSO.pred1[[2]][i,j])) {
    if (is.na(LASSO_pred1[[2]][i,j])) {
      next
    }
    else {
      # LASSO_1[LASSO_res$LASSO.pred1[[2]][i,j],j] <- LASSO_res$LASSO.pred1[[3]][i,j]
      LASSO_1[LASSO_pred1[[2]][i,j],j] <- LASSO_pred1[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:dim(LASSO_pred5[[2]])[1])
  {
    if (is.na(LASSO_pred5[[2]][i,j])) {
      next
    }
    else {
      LASSO_2[LASSO_pred5[[2]][i,j],j] <- LASSO_pred5[[3]][i,j]
    }
  }
  # #####################################################################
  for (i in 1:dim(LASSO_pred2[[2]])[1])
  {
    if (is.na(LASSO_pred2[[2]][i,j])) {
      next
    }
    else {
      LASSO_5[LASSO_pred2[[2]][i,j],j] <- LASSO_pred2[[3]][i,j]
    }
  }
  # #####################################################################
  for (i in 1:dim(LASSO_pred6[[2]])[1])
  {
    if (is.na(LASSO_pred6[[2]][i,j])) {
      next
    }
    else {
      LASSO_6[LASSO_pred6[[2]][i,j],j] <- LASSO_pred6[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:54)
  {
    RF_3[RF_pred3[[2]][i,j],j] <- RF_pred3[[3]][i,j]
    # RF_2[RF_res$RF.pred2[[2]][i,j],j] <- RF_res$RF.pred2[[3]][i,j]
    RF_4[RF_pred7[[2]][i,j],j] <- RF_pred7[[3]][i,j]
  }
  for (i in 1:1122)
  {
    RF_7[RF_pred4[[2]][i,j],j] <- RF_pred4[[3]][i,j]
    RF_8[RF_pred8[[2]][i,j],j] <- RF_pred8[[3]][i,j]
  }
}

#### Calculating rank based on frequency of variable selection for LASSO ####

# Initialize the count holders
counts_FM <- rep(0, dim(LASSO_1)[[1]])
counts_WO <- rep(0, dim(LASSO_2)[[1]])
counts_qFM <- rep(0, dim(LASSO_5)[[1]])
counts_qWO <- rep(0, dim(LASSO_6)[[1]])

# count
for (ii in 1:length(counts_FM)){
  for (j in 1:100){
    if (is.na(LASSO_1[[ii,j]]) || LASSO_1[[ii,j]]==0){
      counts_FM[[ii]] <- counts_FM[[ii]] + 0
    } else {
      counts_FM[[ii]] <- counts_FM[[ii]]+1
    }
    if (is.na(LASSO_2[[ii,j]]) || LASSO_2[[ii,j]]==0){
      counts_WO[[ii]] <- counts_WO[[ii]] + 0
    } else {
      counts_WO[[ii]] <- counts_WO[[ii]]+1
    }
  }
}
for (ii in 1:length(counts_qFM)){
  for (j in 1:100){
    if (is.na(LASSO_5[[ii,j]]) || LASSO_5[[ii,j]]==0){
      counts_qFM[[ii]] <- counts_qFM[[ii]] + 0
    } else {
      counts_qFM[[ii]] <- counts_qFM[[ii]]+1
    }
  }
}
for (ii in 1:length(counts_qWO)){
  for (j in 1:100){
    if (is.na(LASSO_6[[ii,j]]) || LASSO_6[[ii,j]]==0){
      counts_qWO[[ii]] <- counts_qWO[[ii]] + 0
    } else {
      counts_qWO[[ii]] <- counts_qWO[[ii]]+1
    }
  }
}

# convert the count holders to data frams (to have row names)
counts_FM <- data.frame(counts_FM[2:length(counts_FM)], row.names = gverb)
colnames(counts_FM) <- c("Count")
counts_WO <- data.frame(counts_WO[2:length(counts_WO)], row.names = gverb)
colnames(counts_WO) <- c("Count")

counts_qFM <- data.frame(counts_qFM[2:length(counts_qFM)], row.names = colnames(df3))
colnames(counts_qFM) <- c("Count")
counts_qWO <- data.frame(counts_qWO[2:length(counts_qWO)], row.names = colnames(df3))
colnames(counts_qWO) <- c("Count")

# sort and pick out non-zero elements
sorted_FM<- counts_FM[order(-counts_FM$Count), , drop = FALSE]
sorted_WO<- counts_WO[order(-counts_WO$Count), , drop = FALSE]
sorted_qFM<- counts_qFM[order(-counts_qFM$Count), , drop = FALSE]
sorted_qWO<- counts_qWO[order(-counts_qWO$Count), , drop = FALSE]

gg<- match(sorted_FM[,1],0)
cut1 <- which(gg==1)[1]-1
sorted_FM<- sorted_FM[seq(cut1), , FALSE]
gg<- match(sorted_WO[,1],0)
cut2 <- which(gg==1)[1]-1
sorted_WO<- sorted_WO[seq(cut2), , FALSE]

gg<- match(sorted_qFM[,1],0)
cut3 <- which(gg==1)[1]-1
sorted_qFM<- sorted_qFM[seq(cut3), , FALSE]
gg<- match(sorted_qWO[,1],0)
cut4 <- which(gg==1)[1]-1
sorted_qWO<- sorted_qWO[seq(cut4), , FALSE]

# Simple Bar Plot with Added Labels

# set the plot margins to have more space on the left for variable names
op <- par(mar = c(4,20,4,2) + 0.1)

# plot bar chart, rev makes sure the first item is on top (reverse, since default is botton)
barplot(rev(sorted_FM$Count), horiz = TRUE, main="UEFM Variable Frequency", names.arg =rev(rownames(sorted_FM)),las=2)
barplot(rev(sorted_WO$Count), horiz = TRUE, main="Wolf Variable Frequency", names.arg =rev(rownames(sorted_WO)),las=2)

barplot(rev(sorted_qFM$Count), horiz = TRUE, main="Quad UEFM Variable Frequency", names.arg =rev(rownames(sorted_qFM)),las=2)
barplot(rev(sorted_qWO$Count), horiz = TRUE, main="Quad Wolf Variable Frequency", names.arg =rev(rownames(sorted_qWO)),las=2)

# par(par.o) ## reset plot parameters (margins)

## sort LASSO results (use first column for sorting)
# first, find the rows that have at leasst one value
namelist <- rownames(counts_FM)                       # full list of variable names
# g will hold indices of zero frequency variables
g<-rev(order(counts_FM))[1:cut1]
# gnames will contain the names of the "remaining" features (non-zero freq)
gnames <- namelist[c(as.vector(g))]
# remove the first row (intercept) from all three LASSO matrices
LASSO_1 <- data.frame(LASSO_1[-c(1),],row.names = namelist)
# gg will be the data frame without the unwanted features
gg <- data.frame(LASSO_1[c(as.vector(g)),],row.names = gnames)
# mns <- rowMeans(gg, na.rm=TRUE);
FM_Lin_LASSO_sort <- gg#[order(rev(sorted_FM$Count)), , drop = FALSE]

## and Wolf
namelist <- rownames(counts_WO)
g<-rev(order(counts_WO))[1:cut2]
gnames <- namelist[c(as.vector(g))]
LASSO_2 <- data.frame(LASSO_2[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_2[c(as.vector(g)),],row.names = gnames)
Wolf_Lin_LASSO_sort <- gg

## and Quad FM
namelist <- rownames(counts_qFM)
g<-rev(order(counts_qFM))[1:cut3]
gnames <- namelist[c(as.vector(g))]
LASSO_5 <- data.frame(LASSO_5[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_5[c(as.vector(g)),],row.names = gnames)
FM_Quad_LASSO_sort <- gg

## and Quad Wolf
namelist <- rownames(counts_qWO)
g<-rev(order(counts_qWO))[1:cut4]
gnames <- namelist[c(as.vector(g))]
LASSO_6 <- data.frame(LASSO_6[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_6[c(as.vector(g)),],row.names = gnames)
Wolf_Quad_LASSO_sort <- gg

## NOW, sort RF results, easier since no empty cells
FM_Lin_RF_sort <- RF_3[order(-rowMeans(RF_3)), , drop=FALSE]
Wolf_Lin_RF_sort <- RF_4[order(-rowMeans(RF_4)), , drop=FALSE]

FM_Quad_RF_sort <- RF_7[order(-rowMeans(RF_7)), , drop=FALSE]
Wolf_Quad_RF_sort <- RF_8[order(-rowMeans(RF_8)), , drop=FALSE]

#################### plot the results #################################

## Fugl-Meyer:
tit1 <- "UEFM LASSO Coefficients"
tit2 <- "UEFM RF Variable Importance"
tit3 <- "UEFM LASSO Quad Coefs"
tit4 <- "UEFM RF Quad VI"
plot_sideways(FM_Lin_LASSO_sort,0,tit1)
plot_sideways(FM_Lin_RF_sort,25,tit2)
plot_sideways(FM_Quad_LASSO_sort,0,tit3)
plot_sideways(FM_Quad_RF_sort,25,tit4)

## Wolf:
tit1 <- "Wolf LASSO Coefficients"
tit2 <- "Wolf RF Variable Importance"
tit3 <- "Wolf Quad LASSO Coefs"
tit4 <- "Wolf Quad RF VI"

plot_sideways(Wolf_Lin_LASSO_sort,0,tit1)
plot_sideways(Wolf_Lin_RF_sort,25,tit2)
plot_sideways(Wolf_Quad_LASSO_sort,0,tit3)
plot_sideways(Wolf_Quad_RF_sort,25,tit4)


writeMat("R_Lin_results_v6.mat",FM_Lin_LASSO_sort=FM_Lin_LASSO_sort,FM_Lin_RF_sort=FM_Lin_RF_sort,
         Wolf_Lin_LASSO_sort=Wolf_Lin_LASSO_sort,Wolf_Lin_RF_sort=Wolf_Lin_RF_sort,rnam1=rownames(FM_Lin_LASSO_sort),
         rnam2=rownames(FM_Lin_RF_sort),rnam3=rownames(Wolf_Lin_LASSO_sort),
         rnam4=rownames(Wolf_Lin_RF_sort),countFM=sorted_FM$Count, countWO=sorted_WO$Count)
writeMat("countsF_v6.mat",countsFM = counts_FM, countsWO = counts_WO, countsqFM = counts_qFM,
         countsqWO = counts_qWO, countseFM = counts_eFM, countseWO = counts_eWO)
writeMat("Names_v6.mat",
         rnam1=rownames(FM_Lin_LASSO_sort),
         rnam2=rownames(FM_Lin_RF_sort),
         rnam3=rownames(FM_Quad_LASSO_sort),
         rnam4=rownames(FM_Quad_RF_sort),
         rnam5=rownames(Wolf_Lin_LASSO_sort),
         rnam6=rownames(Wolf_Lin_RF_sort),
         rnam7=rownames(Wolf_Quad_LASSO_sort),
         rnam8=rownames(Wolf_Quad_RF_sort),
         rnam9=rownames(FM_NoExtrap_LASSO_sort),
         rnam10=rownames(FM_NoExtrap_RF_sort),
         rnam11=rownames(Wolf_NoExtrap_LASSO_sort),
         rnam12=rownames(Wolf_NoExtrap_RF_sort))
write.csv(FM_Lin_LASSO_sort,"FM_Lin_LASSO_v6.csv")
write.csv(FM_Lin_RF_sort,"FM_Lin_RF_v6.csv")
write.csv(FM_Quad_LASSO_sort,"FM_Quad_LASSO_v6.csv")
write.csv(FM_Quad_RF_sort,"FM_Quad_RF_v6.csv")
write.csv(Wolf_Lin_LASSO_sort,"WO_Lin_LASSO_v6.csv")
write.csv(Wolf_Lin_RF_sort,"WO_Lin_RF_v6.csv")
write.csv(Wolf_Quad_LASSO_sort,"WO_Quad_LASSO_v6.csv")
write.csv(Wolf_Quad_RF_sort,"WO_Quad_RF_v6.csv")

write.csv(FM_NoExtrap_LASSO_sort,"FM_NoExtrap_LASSO_v6.csv")
write.csv(FM_NoExtrap_RF_sort,"FM_NoExtrap_RF_v6.csv")
write.csv(Wolf_NoExtrap_LASSO_sort,"WO_NoExtrap_LASSO_v6.csv")
write.csv(Wolf_NoExtrap_RF_sort,"WO_NoExtrap_RF_v6.csv")

save(FM_Lin_LASSO_sort,FM_Lin_RF_sort,FM_Quad_LASSO_sort,FM_Quad_RF_sort,
     Wolf_Lin_LASSO_sort,Wolf_Lin_RF_sort,Wolf_Quad_LASSO_sort,Wolf_Quad_RF_sort,
     FM_NoExtrap_LASSO_sort,FM_NoExtrap_RF_sort,Wolf_NoExtrap_LASSO_sort,
     Wolf_NoExtrap_RF_sort, counts_eFM, counts_eWO,
     counts_FM,counts_WO,counts_qFM,counts_qWO,file="LS_RF_Results_v6.rda")

#### Discarded Code ####
# 
