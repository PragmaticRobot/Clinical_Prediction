
## This script reads the output from R and then builds bottom-up 
## models to compare the different kinds of models

rm(list = ls()) #clear the workspace
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,e1071,combinat,
               psych, foreach, randomForest, doParallel, inTrees,tableone,caret,mgcv,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)
# options(error=recover)

## set up the parallel pool
registerDoParallel(cores=4) # 4 cores to do the simulations
# load df1 (linear feature set), df3 (quad feature set), yFM, and yWO 
load('2016dec12_Basics.rda')
source("functions.R") # load all functions we might need

## Find the feature order and set up parameters
RF_FM_order <- read.csv('2016dec12_FM_RF_Order.csv')
RF_WO_order <- read.csv('2016dec12_WO_RF_Order.csv')
nSubj   <- 26                  # how many subjects in your dataset?
nFeats  <- 20                  # max number of features on x-axis of final plot
x       <- 1                   # how many do you want to leave out?
xx      <- 1:26                # this will be trashed soon
Folds   <- combn(nSubj,x)      # all the possible combination of Leave x out
nReps   <- dim(Folds)[2]       # make sure to cover all possibilities
nFold   <- 26                  # nFold cross validation for lasso
rm(xx)                         # remove the trash

## sort the top features in decreasing importance, and store the sorted indices
FM_RF_ind <- RF_FM_order$x
WO_RF_ind <- RF_WO_order$x

## Set up placeholders for things
ThisMat_WO_RF     <- as.data.frame(matrix(0,nSubj,nFeats))     # We will add feature to this set

WO_RF_L           <- matrix(0,nReps,nFeats)     # to store the R2 for linear
WO_RF_Q           <- matrix(0,nReps,nFeats)     # to store the R2 for quadratic

ThisBaseQuad      <- matrix(0,nSubj,nFeats)     # Base quadratic set

IsCateg           <- matrix(0,nFeats,1)         # categorical flag (perhaps not needed)
## Remember features 43:51 are categorical
for (i in 1:nFeats)
{
  if (WO_RF_ind[i] >=43 & WO_RF_ind[i] <= 51)
  {
    IsCateg[i] <- 1
  } else
  {
    IsCateg[i] <- 0
  }
}

# Build the models

Out <- yWO;              # what's our output metric here?
for (i in 1:nFeats)
{
  print(i)
  ThisMat_WO_RF[,i] <- as.data.frame(df1[,c(WO_RF_ind[i])])
  colnames(ThisMat_WO_RF)[i] <- colnames(df1[WO_RF_ind[i]])
  for (j in 1:nReps)
  {
    print(j)
    ## For the linear model
    ThisCV <- as.data.frame(ThisMat_WO_RF[,1:i])
    colnames(ThisCV)[1:i] <- colnames(ThisMat_WO_RF)[1:i]
    ThisCV <- ThisCV[-c(Folds[,j]), ,drop = FALSE]
    thisOut <- Out[-c(Folds[,j])] # features 1:i, remove jth subject
    if (i == 1)
    {
      trash = lm(thisOut~data.matrix(ThisCV))
      WO_RF_L[j,i] <- summary(trash)$r.squared
      rm(trash)
    } else
    {
      trash <- cv.glmnet(data.matrix(ThisCV),thisOut,
                         standardize = TRUE, alpha = 1, nlambda = 1000)
      cutoff <- min(trash$cvup)
      means <- trash$cvm
      trash_mod_ind <- max(which(means<cutoff))
      trash2 <- predict(trash, newx=data.matrix(ThisCV),
                        s= trash$lambda[trash_mod_ind])
      trash3 <- lm(thisOut~trash2)
      WO_RF_L[j,i] <- summary(trash3)$r.squared
      rm(trash,trash2,trash3)
    }
    
    ## For the quadratic model
    ThisQuadCV <- SecondOrderThis(ThisCV)
    trash <- cv.glmnet(data.matrix(ThisQuadCV),thisOut, standardize = TRUE,
                       alpha = 1, nlambda = 1000)
    cutoff <- min(trash$cvup)
    means <- trash$cvm
    trash_mod_ind <- max(which(means<cutoff))
    trash2 <- predict(trash, newx=data.matrix(ThisQuadCV),
                      s= trash$lambda[trash_mod_ind])
    trash3 <- lm(thisOut~trash2)
    WO_RF_Q[j,i] <- summary(trash3)$r.squared
    rm(trash,trash2,trash3)
  }
}

save(WO_RF_Q,WO_RF_L,WO_RF_ind,file = '2016dec12_WO_RF_meta_26C1.rda')
write.csv(WO_RF_L,file = '2016dec12_WO_RF_L_R2_meta_26C1.csv')
write.csv(WO_RF_Q,file = '2016dec12_WO_RF_Q_R2_meta_26C1.csv')
