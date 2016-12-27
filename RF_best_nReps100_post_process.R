rm(list = ls())
par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,e1071,
               psych, foreach, randomForest, doParallel, inTrees,tableone,caret,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)

load('2016dec12_Basics.rda')
load('2016dec12_RF_Best_nReps100.rda')

imp <- sort(FM_RF_L_VI[,1], partial=NULL,decreasing=TRUE, index.return=TRUE)

FM_RF_Order <- order(-rowMeans(FM_RF_L_VI))
WO_RF_Order <- order(-rowMeans(WO_RF_L_VI))

rownames(FM_RF_L_VI) <- gverb
rownames(WO_RF_L_VI) <- gverb

FM_Lin_RF_sort <- FM_RF_L_VI[order(-rowMeans(FM_RF_L_VI)), , drop=FALSE]
Wolf_Lin_RF_sort <- WO_RF_L_VI[order(-rowMeans(WO_RF_L_VI)), , drop=FALSE]

write.csv(FM_RF_L_VI, file = '2016dec12_RF_Best_FM_VI.csv')
write.csv(WO_RF_L_VI, file = '2016dec12_RF_Best_WO_VI.csv')

write.csv(FM_RF_Order, file = '2016dec12_FM_RF_Order.csv')
write.csv(WO_RF_Order, file = '2016dec12_WO_RF_Order.csv')

write.csv(FM_Lin_RF_sort, file = '2016dec12_FM_RF_VI_sorted.csv')
write.csv(Wolf_Lin_RF_sort, file = '2016dec12_WO_RF_VI_sorted.csv')
