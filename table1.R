
#########################
## Setup and libraries ##
#########################

# rm(list = ls())
# options(error=recover)
library(MASS)
library(R.matlab)
library(corrplot)
library(GGally)
library(tableone)

#####################################################
## New Wroking folder: Dropbox/Research/NewDec2015 ##
##                 Data Prep - Importing           ##
#####################################################

Data <- readMat("RF_Ready_noscramble.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
# setwd("c:/Users/Yaz/Dropbox/Research/NewDec2015/")
# setup variables
Features <- data.frame(Data$red.DesignMat[1:26,2:56])
Outcome <- data.frame(Data$y.FM,Data$y.PartsFM)
numfeat <- dim(Features)[2]

###############################
## Data Prep - column titles ##
###############################

colnames(Features)[1] <- Data$red.Dnames.full[[2]][[1]][[1,1]]
for (i in 2:numfeat+1)
{
  colnames(Features)[i-1] <- Data$red.Dnames.full[[i]][[1]][[1,1]]
}
colnames(Outcome)[1] <- 'FullFM'
colnames(Outcome)[2] <- 'ArmFM'
CombinedNum <- cbind(Features[,1:42],Features[,53:55],Outcome)
CombinedFull <- cbind(Features,Outcome)
rm(Data,Dpath,i,numfeat)

colnames(Features) <- gsub(" ","",colnames(Features)) # ggpairs doesn't like spaces in varnames
top <- dim(Features)[2]

## We want to use the "Features" Data Frame for the next part(s)

##################################################
## Invert binary variables (there was an error) ##
##################################################
VarNames <- colnames(Features)
varsToFactor <- VarNames[43:52]
Features[varsToFactor] <- lapply(Features[varsToFactor], factor)

for (i in 43:52)
{
  Features[,i] <- as.factor(2-as.integer(Features[,i])) 
}

########################
## Creating Table One ##
########################

tableOne <- CreateTableOne(vars = VarNames, strata = c("EA"), data = Features)

tableOne

