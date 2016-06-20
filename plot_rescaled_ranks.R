## This script plots the fits and calculates RMSE for all the models, this is done for
# rescaled outcome (relative to initial). The function it calls is that same whether it's
# rescaled or not

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

plot_LSRF_fits(InsRVFM,modListFM,rvFM,'Fugl-Meyer')
plot_LSRF_fits(InsRVWO,modListWO,rvWO,'Wolf')

