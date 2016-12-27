
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
LS_FM_order <- read.csv('2016dec12_LS_FM_order.csv')
LS_WO_order <- read.csv('2016dec12_LS_WO_order.csv')
nSubj   <- 26                  # how many subjects in your dataset?
nFeats  <- 20                  # max number of features on x-axis of final plot
x       <- 1                   # how many do you want to leave out?
xx      <- 1:26                # this will be trashed soon
Folds   <- combn(nSubj,x)      # all the possible combination of Leave x out
nReps   <- dim(Folds)[2]       # make sure to cover all possibilities
rm(xx)                         # remove the trash

## sort the top features in decreasing importance, and store the sorted indices
trash <- sort(LS_FM_order$Means, decreasing = TRUE, index.return=TRUE)
FM_LS_ind <- trash$ix
trash <- sort(LS_WO_order$Means, decreasing = TRUE, index.return=TRUE)
WO_LS_ind <- trash$ix
rm(trash)

## Set up placeholders for things
ThisMat_WO_LS     <- as.data.frame(matrix(0,nSubj,nFeats))     # We will add feature to this set
 
WO_LS_L           <- matrix(0,nReps,nFeats)     # to store the R2 for linear
WO_LS_Q           <- matrix(0,nReps,nFeats)     # to store the R2 for quadratic

ThisBaseQuad      <- matrix(0,nSubj,nFeats)     # Base quadratic set

IsCateg           <- matrix(0,nFeats,1)         # categorical flag (perhaps not needed)
## Remember features 43:51 are categorical
for (i in 1:nFeats)
{
  if (WO_LS_ind[i] >=43 & WO_LS_ind[i] <= 51)
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
  ThisMat_WO_LS[,i] <- as.data.frame(df1[,c(WO_LS_ind[i])])
  colnames(ThisMat_WO_LS)[i] <- colnames(df1[WO_LS_ind[i]])
  for (j in 1:nReps)
  {
    print(j)
    ## For the linear model
    ThisCV <- as.data.frame(ThisMat_WO_LS[,1:i])
    colnames(ThisCV)[1:i] <- colnames(ThisMat_WO_LS)[1:i]
    ThisCV <- ThisCV[-c(Folds[,j]), ,drop = FALSE]
    thisOut <- Out[-c(Folds[,j])] # features 1:i, remove jth subject
    if (i == 1)
    {
      trash = lm(thisOut~data.matrix(ThisCV))
      WO_LS_L[j,i] <- summary(trash)$r.squared
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
      WO_LS_L[j,i] <- summary(trash3)$r.squared
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
    WO_LS_Q[j,i] <- summary(trash3)$r.squared
    rm(trash,trash2,trash3)
  }
}

save(WO_LS_Q,WO_LS_L,WO_LS_ind,file = '2016dec12_WO_meta_26C1.rda')
write.csv(WO_LS_L,file = '2016dec12_WO_LS_L_R2_meta_26C1.csv')
write.csv(WO_LS_Q,file = '2016dec12_WO_LS_Q_R2_meta_26C1.csv')
