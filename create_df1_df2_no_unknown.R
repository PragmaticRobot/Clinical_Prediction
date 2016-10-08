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
rm(Data,i,numfeat)

colnames(Features) <- gsub(" ","",colnames(Features)) # ggpairs 
# doesn't like spaces in varnames
top <- dim(Features)[2]

## We want to use the "Features" Data Frame for the next part(s)

############################################################
## Invert binary variables (there was an error in matlab) ##
############################################################
VarNames <- colnames(Features)
varsToFactor <- VarNames[43:52]
# FeaturesF <- Features # a copy of features with binary variables 
# turned into factors (so the original is fine)
Features[varsToFactor] <- lapply(Features[varsToFactor], factor)

for (i in 43:52)
{
  Features[,i] <- as.factor(2-as.integer(Features[,i]))
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

df1 <- df1[,-48]
df2 <- df1
# for (i in 43:52)
# {
#   df2[,i] <- as.integer(df2[,i])-1
#   Features[,i] <- as.integer(Features[,i])-1
#   df1[,i] <- as.integer(df1[,i])-1
# }
df3 <- df1
names2 <- c(colnames(df1))
names3 <- names2

for (i in c(1:dim(df1)[2])) {
  for (j in c(i:dim(df1)[2])){
    if (is.factor(df1[,i]) & is.factor(df1[,j]) & identical(df1[,i],df1[,j])) # If they're the same and both factors, skip
    {
      next
    } else if (is.factor(df1[,i]) & is.factor(df1[,j]) & !identical(df1[,i],df1[,j])) # If not the same and both factors, interact
    {
      name_new <- paste(colnames(df1)[i], ':',colnames(df1)[j], sep="")
      names2 <- c(names2,name_new)
      trash <- interaction(df1[,i],df1[,j],sep = ":")
      df_new <- data.frame(new_name = (trash))
      colnames(df_new) <- name_new
      df2 <- cbind(df2,(df_new))
      # now check for df3, whether any of the "cells" in the table are empty
      if (min(table(df_new))==0){
        next
      } else {
        df3 <- cbind(df3,(df_new))
        names3 <- c(names3,name_new)
      }
    } else if ((is.factor(df1[,i]) & !is.factor(df1[,j])) | (!is.factor(df1[,i]) & is.factor(df1[,j])))
    {
      name_new <- paste(colnames(df1)[i], ':',colnames(df1)[j], sep="")
      names2 <- c(names2,name_new)
      trash <- interaction(df1[,i],df1[,j],sep = ":")
      df_new <- data.frame(new_name = (trash))
      colnames(df_new) <- name_new
      df2 <- cbind(df2,(df_new))
      # now check for df3, whether any of the "cells" in the table are empty
      df3 <- cbind(df3,(df_new))
      names3 <- c(names3,name_new)
    } else
    {
      name_new <- paste(colnames(df1)[i], ':',colnames(df1)[j], sep="")
      names2 <- c(names2,name_new)
      trash <- (df1[,i]*df1[,j])
      df_new <- data.frame(name_new = (trash))
      colnames(df_new) <- name_new
      df2 <- cbind(df2,(df_new))
      # now check for df3, whether any of the "cells" in the table are empty
      df3 <- cbind(df3,(df_new))
      names3 <- c(names3,name_new)
      }
  }
}
colnames(df2) <- c(names2)
colnames(df3) <- c(names3)
# FeaturesQ <- df2               
# # a new data frame containing quadratic predictors with correct column
# # names (thanks Peter!)
# FeaturesQ[varsToFactor] <- lapply(FeaturesQ[varsToFactor], factor)
# 
# for (i in 43:52)
# {
#   FeaturesQ[,i] <- as.factor(2-as.integer(FeaturesQ[,i])) 
# }
