### Let's see what PCA tells us! #########
## This script will pool the speed-related features and run PCA on them,
## see how much they covary

#########################
## Setup and libraries ##
#########################

rm(list = ls()) # clear the environment
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

