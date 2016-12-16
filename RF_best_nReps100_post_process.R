rm(list = ls())
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,e1071,
               psych, foreach, randomForest, doParallel, inTrees,tableone,caret,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)

load('2016dec12_Basics.rda')
load('2016dec12_RF_Best_nReps100.rda')
