
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

Data <- readMat("pqfile.mat", maxLength=NULL, fixNames=TRUE, Verbose=FALSE)
Metrics <- data.frame(Data$AllMets)

x<-NULL
for (i in 1:34){
  x[i]<-Data$FeatNames[[i]][[1]][1]
}
colnames(Metrics) <- x

##################################################
## Invert binary variables (there was an error) ##
##################################################
VarNames <- colnames(Metrics)
varsToFactor <- VarNames[22:31]
Metrics[varsToFactor] <- lapply(Metrics[varsToFactor], factor)

for (i in 22:31)
{
  Metrics[,i] <- as.factor(2-as.integer(Metrics[,i])) 
}

########################
## Creating Table One ##
########################

tableOne <- CreateTableOne(vars = VarNames, strata = c("EA"), data = Metrics)

tableOne

