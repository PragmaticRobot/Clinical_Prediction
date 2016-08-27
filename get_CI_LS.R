# par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,
               psych, foreach, randomForest, doParallel, inTrees,tableone,boot,
               lars,glmnet,coefplot, qpcR, qqman,corrplot,ptw,ggplot2,tcltk2)
source("functions.R")
source("create_df1_df2.R")

registerDoParallel(cores=4) # 4 cores to do the simulations


load("LS_RF_Results.rda")

# Order of variables 1-8
# FM_Lin_LASSO_sort
# FM_Lin_RF_sort
# FM_Quad_LASSO_sort
# FM_Quad_RF_sort
# Wolf_Lin_LASSO_sort
# Wolf_Lin_RF_sort
# Wolf_Quad_LASSO_sort
# Wolf_Quad_RF_sort


Booted1 <- list()
Booted2 <- list()
Booted3 <- list()
Booted4 <- list()
Booted5 <- list()
Booted6 <- list()
Booted7 <- list()
Booted8 <- list()

for (i in 1:17){
  x       = as.double(FM_Lin_LASSO_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted1[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}
for (i in 1:55){
  x       = as.double(FM_Lin_RF_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted2[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}
for (i in 1:33){
  x       = as.double(FM_Quad_LASSO_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted3[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}
for (i in 1:1595){
  x       = as.double(FM_Quad_RF_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted4[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}
for (i in 1:23){
  x       = as.double(Wolf_Lin_LASSO_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted5[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}
for (i in 1:55){
  x       = as.double(Wolf_Lin_RF_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted6[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}
for (i in 1:36){
  x       = as.double(Wolf_Quad_LASSO_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted7[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}
for (i in 1:1595){
  x       = as.double(Wolf_Quad_RF_sort[i,])
  bootmed = apply(matrix(sample(x, rep=TRUE, 10^4*length(x)), nrow=10^4), 1, median, na.rm = TRUE)
  Booted8[[i]] <- quantile(bootmed, c(.025, 0.975), na.rm = TRUE)
}

save(Booted1,Booted2,Booted3,Booted4,Booted5,Booted6,Booted7,Booted8,file="BootedCIs.rda")