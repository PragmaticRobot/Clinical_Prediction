# This version adds bootstrapping to LASSO, so the model predictions
# can be averaged and mean and std can be calculated for these
# predictions, the results will be forwarded to matlab for plotting

#########################
## Setup and libraries ##
#########################

# rm(list = ls()) # clear the environment
library(MASS)
library(R.matlab)

##
LASSO_res <- readMat("lasso_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
RF_res <- readMat("RF_pred.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
load("colnames.rda")
## restructure the data into a data frame
# RF: Random Forests
# LS: LASSO
################################## IMPORTANT NOTE #################################################
######### This section ONLY deals with LINEAR results, quadratic models are not included! #########
###################################################################################################
gm2 <- c("Intercept",gm)    # add "intercept" to the list of feature names for LASSO
gverb <- gm                 # verbose version of feature names

######### Making feature names verbose (looong list) #############
gverb[1] <- "Mean Reaction Time"; gverb[2] <- "Mean Trial Time"
gverb[3] <- "Mean Initial Direction Error"; gverb[4] <- "Mean Pre-Movement Speed"
gverb[5] <- "Mean Max Speed"; gverb[6] <- "Mean Hand Path Length"
gverb[7] <- "Mean Initial Movement Ratio"; gverb[8] <- "Mean Speed Ratio"
gverb[9] <- "Mean Path Length Ratio"; gverb[10] <- "Mean Number of Speed Peaks"
gverb[11] <- "Mean Max Perpendicular Distance"; gverb[12] <- "Mean Percent of Movement in the Target Direction"
gverb[13] <- "Mean Mean Arrest Period Ratio"; gverb[14] <- "Max Reaction Time"
gverb[15] <- "Max Trial Time"; gverb[16] <- "Max Initial Direction Error"
gverb[17] <- "Max Pre-Movement Speed"; gverb[18] <- "Max Max Speed"
gverb[19] <- "Max Hand Path Length"; gverb[20] <- "Max Initial Movement Ratio"
gverb[21] <- "Max Path Length Ratio"; gverb[22] <- "Max Number of Speed Peaks"
gverb[23] <- "Max Max Perpendicular Distance"; gverb[24] <- "Max Percent of Movement in the Target Direction"
gverb[25] <- "Max Mean Arrest Period Ratio"; gverb[26] <- "Var Reaction Time"
gverb[27] <- "Var Trial Time"; gverb[28] <- "Var Initial Direction Error"
gverb[29] <- "Var Pre-Movement Speed"; gverb[30] <- "Var Max Speed"
gverb[31] <- "Var Hand Path Length"; gverb[32] <- "Var Initial Movement Ratio"
gverb[33] <- "Var Speed Ratio"; gverb[34] <- "Var Path Length Ratio"
gverb[35] <- "Var Number of Speed Peaks"; gverb[36] <- "Var Max Perpendicular Distance"
gverb[37] <- "Var Percent of Movement in the Target Direction"; gverb[38] <- "Var Mean Arrest Period Ratio"
gverb[39] <- "Age"; gverb[40] <- "Height"; gverb[41] <- "Mass"; 
gverb[42] <- "Months Post-Stroke"; gverb[43] <- "Sex"; gverb[44] <- "Left Hand Dominant?"
gverb[45] <- "Left Side Affected?"; gverb[46] <- "Dominant Side Affected?"
gverb[47] <- "Hemorrhagic Stroke?"; gverb[48] <- "Unknown Stroke Location?"
gverb[49] <- "Cortical Stroke?"; gverb[50] <- "Subcortical Stroke?"
gverb[51] <- "Brainstem Stroke?"; gverb[52] <- "Error-Augmentation Treatment?"
gverb[53] <- "Initial Fugl-Meyer Score"; gverb[54] <- "Initial Wolf Motor Function Score"
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
counts_FM <- data.frame(counts_FM[2:length(counts_FM)], row.names = gverb)
colnames(counts_FM) <- c("Count")
counts_pFM <- data.frame(counts_pFM[2:length(counts_pFM)], row.names = gverb)
colnames(counts_pFM) <- c("Count")
counts_WO <- data.frame(counts_WO[2:length(counts_WO)], row.names = gverb)
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
png("UEFM_counts.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
barplot(rev(sorted_FM$Count), horiz = TRUE, main="UEFM Variable Frequency", names.arg =rev(rownames(sorted_FM)),las=2)
dev.off()

png("ArmFM_counts.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
barplot(rev(sorted_pFM$Count), horiz = TRUE, main="Arm-Only FM Variable Frequency", names.arg =rev(rownames(sorted_pFM)),las=2)
dev.off()

png("Wolf_counts.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
barplot(rev(sorted_WO$Count), horiz = TRUE, main="Wolf Variable Frequency", names.arg =rev(rownames(sorted_WO)),las=2)
dev.off()

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

png("UEFM_LASSO.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(FM_Lin_LASSO_sort,0,tit1)
dev.off()

png("UEFM_RF.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(FM_Lin_RF_sort,15,tit2)
dev.off()

## Arm-Only FM:
tit1 <- "Arm-Only FM LASSO Coefficients"
tit2 <- "Arm-Only FM RF Variable Importance"
png("ArmFM_LASSO.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(ArmFM_Lin_LASSO_sort,0,tit1)
dev.off()

png("ArmFM_RF.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(ArmFM_Lin_RF_sort,15,tit2)
dev.off()

## Wolf:
tit1 <- "Wolf LASSO Coefficients"
tit2 <- "Wolf RF Variable Importance"
png("Wolf_LASSO.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(Wolf_Lin_LASSO_sort,0,tit1)
dev.off()

png("Wolf_RF.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
plot_sideways(Wolf_Lin_RF_sort,15,tit2)
dev.off()

########################### Functions ######################################
#
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
  offs <- runif(length(x),-0.15,0.15)
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
    offs <- runif(length(x),-0.15,0.15)
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


##################################