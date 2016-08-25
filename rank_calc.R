# This version adds bootstrapping to LASSO, so the model predictions
# can be averaged and mean and std can be calculated for these
# predictions, the results will be forwarded to matlab for plotting

# rm(list = ls()) # clear the environment

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

####################################
LASSO_res <- readMat("lasso_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
RF_res <- readMat("RF_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
load("colnames.rda")
## restructure the data into a data frame
# RF: Random Forests
# LS: LASSO
# initialize data frames with the correct sizes and add row names for features
colnmquad <- c("Intercept",colnames(df2)) #colnames for quadratic models

LASSO_1 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
# LASSO_2 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
LASSO_3 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
LASSO_4 <- data.frame(matrix(vector(),1596,100),row.names = colnmquad)
LASSO_6 <- data.frame(matrix(vector(),1596,100),row.names = colnmquad)
RF_1 <- data.frame(matrix(vector(),55,100), row.names = gverb)
# RF_2 <- data.frame(matrix(vector(),55,100), row.names = gverb)
RF_3 <- data.frame(matrix(vector(),55,100), row.names = gverb)
RF_4 <- data.frame(matrix(vector(),1595,100), row.names = colnames(df2))
RF_6 <- data.frame(matrix(vector(),1595,100), row.names = colnames(df2))

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
  # for (i in 1:dim(LASSO_res$LASSO.pred2[[2]])[1])
  # {
  #   if (is.na(LASSO_res$LASSO.pred2[[2]][i,j])) {
  #     next 
  #   }
  #   else {
  #     LASSO_2[LASSO_res$LASSO.pred2[[2]][i,j],j] <- LASSO_res$LASSO.pred2[[3]][i,j]
  #   }
  # }
  # #####################################################################
  for (i in 1:dim(LASSO_res$LASSO.pred3[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred3[[2]][i,j])) {
      next 
    }
    else {
      LASSO_3[LASSO_res$LASSO.pred3[[2]][i,j],j] <- LASSO_res$LASSO.pred3[[3]][i,j]
    }
  }
  # #####################################################################
  for (i in 1:dim(LASSO_res$LASSO.pred4[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred4[[2]][i,j])) {
      next 
    }
    else {
      LASSO_4[LASSO_res$LASSO.pred4[[2]][i,j],j] <- LASSO_res$LASSO.pred4[[3]][i,j]
    }
  }
  # #####################################################################
  for (i in 1:dim(LASSO_res$LASSO.pred6[[2]])[1])
  {
    if (is.na(LASSO_res$LASSO.pred6[[2]][i,j])) {
      next 
    }
    else {
      LASSO_6[LASSO_res$LASSO.pred6[[2]][i,j],j] <- LASSO_res$LASSO.pred6[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:55)
  {
    RF_1[RF_res$RF.pred1[[2]][i,j],j] <- RF_res$RF.pred1[[3]][i,j]
    # RF_2[RF_res$RF.pred2[[2]][i,j],j] <- RF_res$RF.pred2[[3]][i,j]
    RF_3[RF_res$RF.pred3[[2]][i,j],j] <- RF_res$RF.pred3[[3]][i,j]
  }
  for (i in 1:1595)
  {
    RF_4[RF_res$RF.pred4[[2]][i,j],j] <- RF_res$RF.pred4[[3]][i,j]
    # RF_2[RF_res$RF.pred2[[2]][i,j],j] <- RF_res$RF.pred2[[3]][i,j]
    RF_6[RF_res$RF.pred6[[2]][i,j],j] <- RF_res$RF.pred6[[3]][i,j]
  }
}

## Calculating rank based on frequency of variable selection for LASSO

# Initialize the count holders
counts_FM <- rep(0, dim(LASSO_1)[[1]])
# counts_pFM <- rep(0, dim(LASSO_2)[[1]])
counts_WO <- rep(0, dim(LASSO_3)[[1]])
counts_qFM <- rep(0, dim(LASSO_4)[[1]])
counts_qWO <- rep(0, dim(LASSO_6)[[1]])

# count
for (ii in 1:length(counts_FM)){
  for (j in 1:100){
    if (is.na(LASSO_1[[ii,j]]) || LASSO_1[[ii,j]]==0){
      counts_FM[[ii]] <- counts_FM[[ii]] + 0
    } else {
      counts_FM[[ii]] <- counts_FM[[ii]]+1
    }
    # if (is.na(LASSO_2[[ii,j]]) || LASSO_2[[ii,j]]==0){
    #   counts_pFM[[ii]] <- counts_pFM[[ii]] + 0
    # } else {
    #   counts_pFM[[ii]] <- counts_pFM[[ii]]+1
    # }
    if (is.na(LASSO_3[[ii,j]]) || LASSO_3[[ii,j]]==0){
      counts_WO[[ii]] <- counts_WO[[ii]] + 0
    } else {
      counts_WO[[ii]] <- counts_WO[[ii]]+1
    }
  }
}
for (ii in 1:1596){
  for (j in 1:100){
    if (is.na(LASSO_4[[ii,j]]) || LASSO_4[[ii,j]]==0){
      counts_qFM[[ii]] <- counts_qFM[[ii]] + 0
    } else {
      counts_qFM[[ii]] <- counts_qFM[[ii]]+1
    }
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
# counts_pFM <- data.frame(counts_pFM[2:length(counts_pFM)], row.names = gverb)
# colnames(counts_pFM) <- c("Count")
counts_WO <- data.frame(counts_WO[2:length(counts_WO)], row.names = gverb)
colnames(counts_WO) <- c("Count")

counts_qFM <- data.frame(counts_qFM[2:length(counts_qFM)], row.names = colnames(df2))
colnames(counts_qFM) <- c("Count")
counts_qWO <- data.frame(counts_qWO[2:length(counts_qWO)], row.names = colnames(df2))
colnames(counts_qWO) <- c("Count")

# sort and pick out non-zero elements
sorted_FM<- counts_FM[order(-counts_FM$Count), , drop = FALSE]
# sorted_pFM<- counts_pFM[order(-counts_pFM$Count), , drop = FALSE]
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
# barplot(rev(sorted_pFM$Count), horiz = TRUE, main="Arm-Only FM Variable Frequency", names.arg =rev(rownames(sorted_pFM)),las=2)
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
LASSO_3 <- data.frame(LASSO_3[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_3[c(as.vector(g)),],row.names = gnames)
Wolf_Lin_LASSO_sort <- gg

## and Quad FM
namelist <- rownames(counts_qFM)
g<-rev(order(counts_qFM))[1:cut3]
gnames <- namelist[c(as.vector(g))]
LASSO_4 <- data.frame(LASSO_4[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_4[c(as.vector(g)),],row.names = gnames)
FM_Quad_LASSO_sort <- gg

## and Quad Wolf
namelist <- rownames(counts_qWO)
g<-rev(order(counts_qWO))[1:cut4]
gnames <- namelist[c(as.vector(g))]
LASSO_6 <- data.frame(LASSO_6[-c(1),],row.names = namelist)
gg <- data.frame(LASSO_6[c(as.vector(g)),],row.names = gnames)
Wolf_Quad_LASSO_sort <- gg

# ## now repeat for Arm FM
# g<-which(counts_pFM==0,arr.ind = T)
# gnames <- namelist[-c(as.vector(g[,1]))]
# LASSO_2 <- data.frame(LASSO_2[-c(1),],row.names = namelist)
# gg <- data.frame(LASSO_2[-c(as.vector(g[,1])),],row.names = gnames)
# mns <- rowMeans(gg, na.rm=TRUE);
# ArmFM_Lin_LASSO_sort <- gg[order(-abs(mns)), , drop = FALSE]

## NOW, sort RF results, easier since no empty cells
FM_Lin_RF_sort <- RF_1[order(-RF_1$X1), , drop=FALSE]
# ArmFM_Lin_RF_sort <- RF_2[order(-RF_2$X1), , drop=FALSE]
Wolf_Lin_RF_sort <- RF_3[order(-RF_3$X1), , drop=FALSE]
FM_Quad_RF_sort <- RF_4[order(-RF_4$X1), , drop=FALSE]
Wolf_Quad_RF_sort <- RF_6[order(-RF_6$X1), , drop=FALSE]

#################### plot the results #################################

## Fugl-Meyer:
tit1 <- "UEFM LASSO Coefficients"
tit2 <- "UEFM RF Variable Importance"
tit3 <- "UEFM LASSO Quad Coefs"
tit4 <- "UEFM RF Quad VI"
plot_sideways(FM_Lin_LASSO_sort,0,tit1)
plot_sideways(FM_Lin_RF_sort,15,tit2)
plot_sideways(FM_Quad_LASSO_sort,0,tit3)
plot_sideways(FM_Quad_RF_sort,15,tit4)

# ## Arm-Only FM:
# tit1 <- "Arm-Only FM LASSO Coefficients"
# tit2 <- "Arm-Only FM RF Variable Importance"
# # png("ArmFM_LASSO.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
# plot_sideways(ArmFM_Lin_LASSO_sort,0,tit1)
# # dev.off()
# plot_sideways(ArmFM_Lin_RF_sort,15,tit2)

## Wolf:
tit1 <- "Wolf LASSO Coefficients"
tit2 <- "Wolf RF Variable Importance"
tit3 <- "Wolf Quad LASSO Coefs"
tit4 <- "Wolf Quad RF VI"
plot_sideways(Wolf_Lin_LASSO_sort,0,tit1)
plot_sideways(Wolf_Lin_RF_sort,15,tit2)
plot_sideways(Wolf_Quad_LASSO_sort,0,tit1)
plot_sideways(Wolf_Quad_RF_sort,15,tit2)

# writeMat("R_results.mat",FM_Lin_LASSO_sort=FM_Lin_LASSO_sort,FM_Lin_RF_sort=FM_Lin_RF_sort,
#          Wolf_Lin_LASSO_sort=Wolf_Lin_LASSO_sort,Wolf_Lin_RF_sort=Wolf_Lin_RF_sort,rnam1=rownames(FM_Lin_LASSO_sort),
#          rnam2=rownames(FM_Lin_RF_sort),rnam3=rownames(Wolf_Lin_LASSO_sort),
#          rnam4=rownames(Wolf_Lin_RF_sort),countFM=sorted_FM$Count, countWO=sorted_WO$Count)
writeMat("countsF.mat",countsFM = counts_FM, countsWO = counts_WO, countsqFM = counts_qFM,
         countsqWO = counts_qWO)
writeMat("Names.mat",
         rnam1=rownames(FM_Lin_LASSO_sort),
         rnam2=rownames(FM_Lin_RF_sort),
         rnam3=rownames(FM_Quad_LASSO_sort),
         rnam4=rownames(FM_Quad_RF_sort),
         rnam5=rownames(Wolf_Lin_LASSO_sort),
         rnam6=rownames(Wolf_Lin_RF_sort),
         rnam7=rownames(Wolf_Quad_LASSO_sort),
         rnam8=rownames(Wolf_Quad_RF_sort))
write.csv(FM_Lin_LASSO_sort,"FM_Lin_LASSO.csv")
write.csv(FM_Lin_RF_sort,"FM_Lin_RF.csv")
write.csv(FM_Quad_LASSO_sort,"FM_Quad_LASSO.csv")
write.csv(FM_Quad_RF_sort,"FM_Quad_RF.csv")
write.csv(Wolf_Lin_LASSO_sort,"WO_Lin_LASSO.csv")
write.csv(Wolf_Lin_RF_sort,"WO_Lin_RF.csv")
write.csv(Wolf_Quad_LASSO_sort,"WO_Quad_LASSO.csv")
write.csv(Wolf_Quad_RF_sort,"WO_Quad_RF.csv")

save(FM_Lin_LASSO_sort,FM_Lin_RF_sort,FM_Quad_LASSO_sort,FM_Quad_RF_sort,
     Wolf_Lin_LASSO_sort,Wolf_Lin_RF_sort,Wolf_Quad_LASSO_sort,Wolf_Quad_RF_sort,
     counts_FM,counts_WO,counts_qFM,counts_qWO,file="LS_RF_Results.rda")