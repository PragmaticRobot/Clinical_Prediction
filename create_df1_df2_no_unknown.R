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

# gm2 <- c("Intercept",gm)    # add "intercept" to the list of feature names for LASSO
# gverb <- gm                 # verbose version of feature names
gverb <- NULL
######### Making feature names verbose (looong list) #############
gverb[1] <- "Mean Reaction Time"; gverb[2] <- "Mean Trial Time"
gverb[3] <- "Mean Initial Direction Error"; gverb[4] <- "Mean Pre-Movement Speed"
gverb[5] <- "Mean Max Speed"; # gverb[6] <- "Mean Hand Path Length"
gverb[6] <- "Mean Initial Movement Ratio"; gverb[7] <- "Mean Speed Ratio"
gverb[8] <- "Mean Path Length Ratio"; gverb[9] <- "Mean Number of Speed Peaks"
gverb[10] <- "Mean Max Perpendicular Distance"; 
gverb[11] <- "Mean Percent of Movement in the Target Direction"
gverb[12] <- "Mean Mean Arrest Period Ratio"; gverb[13] <- "Max Reaction Time"
gverb[14] <- "Max Trial Time"; gverb[15] <- "Max Initial Direction Error"
gverb[16] <- "Max Pre-Movement Speed"; gverb[17] <- "Max Max Speed"
# gverb[18] <- "Max Hand Path Length"; 
gverb[18] <- "Max Initial Movement Ratio"
gverb[19] <- "Max Path Length Ratio"; gverb[20] <- "Max Number of Speed Peaks"
gverb[21] <- "Max Max Perpendicular Distance"; 
gverb[22] <- "Max Percent of Movement in the Target Direction"
gverb[23] <- "Max Mean Arrest Period Ratio"; gverb[24] <- "Var Reaction Time"
gverb[25] <- "Var Trial Time"; gverb[26] <- "Var Initial Direction Error"
gverb[27] <- "Var Pre-Movement Speed"; gverb[28] <- "Var Max Speed"
# gverb[31] <- "Var Hand Path Length"; 
gverb[29] <- "Var Initial Movement Ratio"
gverb[30] <- "Var Speed Ratio"; gverb[31] <- "Var Path Length Ratio"
gverb[32] <- "Var Number of Speed Peaks"; 
gverb[33] <- "Var Max Perpendicular Distance"
gverb[34] <- "Var Percent of Movement in the Target Direction"; 
gverb[35] <- "Var Mean Arrest Period Ratio"
gverb[36] <- "Age"; gverb[37] <- "Height"; gverb[38] <- "Mass"; 
gverb[39] <- "Months Post-Stroke"; gverb[40] <- "Sex"; 
gverb[41] <- "Left Hand Dominant?"
gverb[42] <- "Left Side Affected?"; gverb[43] <- "Dominant Side Affected?"
gverb[44] <- "Hemorrhagic Stroke?"; # gverb[45] <- "Unknown Stroke Location?"
gverb[45] <- "Cortical Stroke?"; gverb[46] <- "Subcortical Stroke?"
gverb[47] <- "Brainstem Stroke?"; gverb[48] <- "Error-Augmentation Treatment?"
gverb[49] <- "Initial Fugl-Meyer Score"; 
gverb[50] <- "Initial Wolf Motor Function Score"
gverb[51] <- "Initial Box-and-Blocks Score"

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

df1 <- df1[,-48] # remove (unknown stroke location column)

# to remove all hand path length features (3 columns), uncomment next line
# columns to remove: 6, 19, 31, removed one at a time would be 6, 18, 29
df1 <- df1[, -6] # remove mean HPL
df1 <- df1[, -18] # remove max HPL
df1 <- df1[, -29] # remove var HPL

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
save(df1,df3,yFM,yWO,gverb,gverb2, file = '2016dec12_Basics.rda')