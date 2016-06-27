## This script plots the fits and calculates RMSE for all the models, this is done for
# rescaled outcome (relative to initial). The function it calls is that same whether it's
# rescaled or not

<<<<<<< HEAD
# rm(list = ls()) # clear the environment
=======
rm(list = ls()) # clear the environment
>>>>>>> 5b757334244bc6316360ec2f75f9d47f18b5f4a2

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

load("rv_cv_models.rda")
load("rescaled_out.rda")
load("rv_cv_all_LSRF.rda")

# build the linear models
rvfm_lin_lasso_mod <- lm(rowMeans(rvFM_lin_lasso)~rvFM)
rvfm_quad_lasso_mod <- lm(rowMeans(rvFM_quad_lasso)~rvFM)
rvfm_lin_rf_mod <- lm(rowMeans(rvFM_lin_RF)~rvFM)
rvfm_quad_rf_mod <- lm(rowMeans(rvFM_quad_RF)~rvFM)

rvwo_lin_lasso_mod <- lm(rowMeans(rvWO_lin_lasso)~rvWO)
rvwo_quad_lasso_mod <- lm(rowMeans(rvWO_quad_lasso)~rvWO)
rvwo_lin_rf_mod <- lm(rowMeans(rvWO_lin_RF)~rvWO)
rvwo_quad_rf_mod <- lm(rowMeans(rvWO_quad_RF)~rvWO)

# put all the inputs (rowMeans) in a matrix for passing to the plotting function

InsRVFM <- cbind(rowMeans(rvFM_lin_lasso),rowMeans(rvFM_quad_lasso),
                 rowMeans(rvFM_lin_RF), rowMeans(rvFM_quad_RF)) # this is for rvFM

InsRVWO <- cbind(rowMeans(rvWO_lin_lasso),rowMeans(rvWO_quad_lasso),
                 rowMeans(rvWO_lin_RF), rowMeans(rvWO_quad_RF)) # this is for rvFM

modListFM <- list(rvfm_lin_lasso_mod,rvfm_quad_lasso_mod,
                  rvfm_lin_rf_mod,rvfm_quad_rf_mod)
modListWO <- list(rvwo_lin_lasso_mod,rvwo_quad_lasso_mod,
                  rvwo_lin_rf_mod,rvwo_quad_rf_mod)

plot_LSRF_fits(InsRVFM,modListFM,1,rvFM,'Fugl-Meyer')
plot_LSRF_fits(InsRVWO,modListWO,1,rvWO,'Wolf')

# next we arrange things to calculate ranks. 1 is FM, 2 is WO
LASSO_1 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
LASSO_2 <- data.frame(matrix(vector(),56,100), row.names = gverb2)
RF_1 <- data.frame(matrix(vector(),55,100), row.names = gverb)
RF_2 <- data.frame(matrix(vector(),55,100), row.names = gverb)

for (j in 1:100) # loop over cross-validations
{
  for (i in 1:dim(LASSO_pred_rv1[[2]])[1])
  {
    if (is.na(LASSO_pred_rv1[[2]][i,j])) {
      next
    }
    else {
      LASSO_1[LASSO_pred_rv1[[2]][i,j],j] <- LASSO_pred_rv1[[3]][i,j]
    }
  }
  ####################################################################
  for (i in 1:dim(LASSO_pred_rv2[[2]])[1])
  {
    if (is.na(LASSO_pred_rv2[[2]][i,j])) {
      next
    }
    else {
      LASSO_2[LASSO_pred_rv2[[2]][i,j],j] <- LASSO_pred_rv2[[3]][i,j]
    }
  }
  #####################################################################
  for (i in 1:55)
  {
    RF_1[RF_pred_rv1[[2]][i,j],j] <- RF_pred_rv1[[3]][i,j]
    RF_2[RF_pred_rv2[[2]][i,j],j] <- RF_pred_rv2[[3]][i,j]
  }
}

## Calculating rank based on frequency of variable selection for LASSO

# Initialize the count holders
counts_FM <- rep(0, dim(LASSO_1)[[1]])
counts_WO <- rep(0, dim(LASSO_2)[[1]])

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

# convert the count holders to data frams (to have row names)
counts_FM <- data.frame(counts_FM[2:length(counts_FM)], row.names = gverb)
colnames(counts_FM) <- c("Count")
counts_WO <- data.frame(counts_WO[2:length(counts_WO)], row.names = gverb)
colnames(counts_WO) <- c("Count")

# sort and pick out non-zero elements
sorted_FM<- counts_FM[order(-counts_FM$Count), , drop = FALSE]
sorted_WO<- counts_WO[order(-counts_WO$Count), , drop = FALSE]
gg<- match(sorted_FM[,1],0)
cut <- which(gg==1)[1]-1
sorted_FM<- sorted_FM[seq(cut), , FALSE]

## In this case, the LASSO wolf models failed completely. 
# next steps will only be done with the non-failed models
# Linear LASSO FM
# Linear RF FM
# Linear RF WO

# Simple Bar Plot with Added Labels

# set the plot margins to have more space on the left for variable names
op <- par(mar = c(4,20,4,2) + 0.1)

# plot bar chart, rev makes sure the first item is on top (reverse, since default is botton)
# png("UEFM_counts.png",width=1500, height = 1014, units = "px", pointsize = 11, bg="white", res=72)
barplot(rev(sorted_FM$Count), horiz = TRUE, main="UEFM Variable Frequency", names.arg =rev(rownames(sorted_FM)),las=2)
# dev.off()

par(par.o) ## reset plot parameters (margins)

## sort LASSO results (use first column for sorting)
# first, find the rows that have at leasst one value
namelist <- rownames(counts_FM)                       # full list of variable names
# g will hold indices of zero frequency variables
g<-rev(order(counts_FM))[1:cut]
# gnames will contain the names of the "remaining" features (non-zero freq)
gnames <- namelist[c(as.vector(g))]
# remove the first row (intercept) from all three LASSO matrices
LASSO_1 <- data.frame(LASSO_1[-c(1),],row.names = namelist)
# gg will be the data frame without the unwanted features
gg <- data.frame(LASSO_1[c(as.vector(g)),],row.names = gnames)
# mns <- rowMeans(gg, na.rm=TRUE);
FM_Lin_LASSO_sort <- gg#[order(rev(sorted_FM$Count)), , drop = FALSE]

## NOW, sort RF results, easier since no empty cells
FM_Lin_RF_sort <- RF_1[order(-RF_1$X1), , drop=FALSE]
Wolf_Lin_RF_sort <- RF_2[order(-RF_2$X1), , drop=FALSE]

#################### plot the results #################################

## Fugl-Meyer:
tit1 <- "UEFM LASSO Coefficients"
tit2 <- "UEFM RF Variable Importance"
plot_sideways(FM_Lin_LASSO_sort,0,tit1)
plot_sideways(FM_Lin_RF_sort,15,tit2)

## Wolf:
tit2 <- "Wolf RF Variable Importance"
plot_sideways(Wolf_Lin_RF_sort,15,tit2)
