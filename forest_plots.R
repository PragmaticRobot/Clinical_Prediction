#### Prep ####

# par.o <- par()
library(pacman)
pacman::p_load(MASS, R.matlab, devtools,rgl, nFactors, FactoMineR,GGally,forestplot,
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

##### Initialize empty vectors to hold data (LowCI, highCI) ####
low1 <- vector()
low2 <- vector()
low3 <- vector()
low4 <- vector()
low5 <- vector()
low6 <- vector()
low7 <- vector()
low8 <- vector()

high1 <- vector()
high2 <- vector()
high3 <- vector()
high4 <- vector()
high5 <- vector()
high6 <- vector()
high7 <- vector()
high8 <- vector()

#### Populate the low/high vectors ####
for (i in 1:length(Booted1)){
  low1[i] <- as.double(Booted1[[i]][1])
  high1[i] <- as.double(Booted1[[i]][2])
}

for (i in 1:length(Booted2)){
  low2[i] <- as.double(Booted2[[i]][1])
  high2[i] <- as.double(Booted2[[i]][2])
}

for (i in 1:length(Booted3)){
  low3[i] <- as.double(Booted3[[i]][1])
  high3[i] <- as.double(Booted3[[i]][2])
}

for (i in 1:length(Booted4)){
  low4[i] <- as.double(Booted4[[i]][1])
  high4[i] <- as.double(Booted4[[i]][2])
}

for (i in 1:length(Booted5)){
  low5[i] <- as.double(Booted5[[i]][1])
  high5[i] <- as.double(Booted5[[i]][2])
}

for (i in 1:length(Booted6)){
  low6[i] <- as.double(Booted6[[i]][1])
  high6[i] <- as.double(Booted6[[i]][2])
}

for (i in 1:length(Booted7)){
  low7[i] <- as.double(Booted7[[i]][1])
  high7[i] <- as.double(Booted7[[i]][2])
}

for (i in 1:length(Booted8)){
  low8[i] <- as.double(Booted8[[i]][1])
  high8[i] <- as.double(Booted8[[i]][2])
}

#### Extract variable names from LS_RF_Results for plotting ####
rnam1 <- rownames(FM_Lin_LASSO_sort)
rnam2 <- rownames(FM_Lin_RF_sort)
rnam3 <- rownames(FM_Quad_LASSO_sort)
rnam4 <- rownames(FM_Quad_RF_sort)
rnam5 <- rownames(Wolf_Lin_LASSO_sort)
rnam6 <- rownames(Wolf_Lin_RF_sort)
rnam7 <- rownames(Wolf_Quad_LASSO_sort)
rnam8 <- rownames(Wolf_Quad_RF_sort)

#### Get means and medians of feature sets ####

colnames(df1) <- rownames(counts_FM) 

# mnsdf1 <- as.data.frame(apply(as.matrix(df1), 2, mean))
# colnames(mnsdf1) <- "means"
# mnsdf2 <- as.data.frame(apply(as.matrix(df2), 2, mean))
# colnames(mnsdf2) <- "means"

#### Notes at this point: ####
# df1 contains original feature set
# df2 contains quadratic feature set
# we have means, medians, lowCIs, highCIs
# we use the means for scaling purposes
# forest plot the medians and CIs

#### Forest plot ####

# Order of variables 1-8
# FM_Lin_LASSO_sort
# FM_Lin_RF_sort
# FM_Quad_LASSO_sort
# FM_Quad_RF_sort
# Wolf_Lin_LASSO_sort
# Wolf_Lin_RF_sort
# Wolf_Quad_LASSO_sort
# Wolf_Quad_RF_sort

mids1 <- (low1 + high1)/2
mids2 <- (low2 + high2)/2
mids3 <- (low3 + high3)/2
mids4 <- (low4 + high4)/2
mids5 <- (low5 + high5)/2
mids6 <- (low6 + high6)/2
mids7 <- (low7 + high7)/2
mids8 <- (low8 + high8)/2

#### For Fugl-Meyer ####
tabletext<-cbind(
  c("Feature",rnam1[1:15]),
  c(rep(NA,16)))

forestplot(tabletext,
           mean = c(NA,mids1[1:15]),
           lower = c(NA, low1[1:15]), 
           upper = c(NA,high1[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-9,14),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Coefficient",lwd.ci = 1.4,boxsize = 0.1)

tabletext<-cbind(
  c("Feature",rnam2[1:15]),
  c(rep(NA,16)))
max2 <- max(FM_Lin_RF_sort)
low2 <- low2/max2
high2 <- high2/max2
mids2 <- mids2/max2
forestplot(tabletext,
           mean = c(NA,mids2[1:15]),
           lower = c(NA, low2[1:15]), 
           upper = c(NA,high2[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-0.1,1),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Relative Feature importance",lwd.ci = 1.4,boxsize = 0.1)

tabletext<-cbind(
  c("Feature",rnam3[1:15]),
  c(rep(NA,16)))
forestplot(tabletext,
           mean = c(NA,mids3[1:15]),
           lower = c(NA, low3[1:15]), 
           upper = c(NA,high3[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-16,26),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Coefficient",lwd.ci = 1.4,boxsize = 0.1)

tabletext<-cbind(
  c("Feature",rnam4[1:15]),
  c(rep(NA,16)))
max4 <- max(FM_Quad_RF_sort)
low4 <- low4/max4
high4 <- high4/max4
mids4 <- mids4/max4
forestplot(tabletext,
           mean = c(NA,mids4[1:15]),
           lower = c(NA, low4[1:15]), 
           upper = c(NA,high4[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-0.1,1),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Relative Feature importance",lwd.ci = 1.4,boxsize = 0.1)

#### For Wolf ####

tabletext<-cbind(
  c("Feature",rnam5[1:15]),
  c(rep(NA,16)))

forestplot(tabletext,
           mean = c(NA,mids5[1:15]),
           lower = c(NA, low5[1:15]), 
           upper = c(NA,high5[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-9,14),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Coefficient",lwd.ci = 1.4,boxsize = 0.1)

tabletext<-cbind(
  c("Feature",rnam6[1:15]),
  c(rep(NA,16)))
max6 <- max(Wolf_Lin_RF_sort)
low6 <- low6/max6
high6 <- high6/max6
mids6 <- mids6/max6
forestplot(tabletext,
           mean = c(NA,mids6[1:15]),
           lower = c(NA, low6[1:15]), 
           upper = c(NA,high6[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-0.1,1),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Relative Feature importance",lwd.ci = 1.4,boxsize = 0.1)

tabletext<-cbind(
  c("Feature",rnam7[1:15]),
  c(rep(NA,16)))
forestplot(tabletext,
           mean = c(NA,mids7[1:15]),
           lower = c(NA, low7[1:15]), 
           upper = c(NA,high7[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-16,26),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Coefficient",lwd.ci = 1.4,boxsize = 0.1)

tabletext<-cbind(
  c("Feature",rnam8[1:15]),
  c(rep(NA,16)))
max8 <- max(Wolf_Quad_RF_sort)
low8 <- low8/max8
high8 <- high8/max8
mids8 <- mids8/max8
forestplot(tabletext,
           mean = c(NA,mids8[1:15]),
           lower = c(NA, low8[1:15]), 
           upper = c(NA,high8[1:15]),
           new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,16)),
           clip=c(-0.1,1),
           xlog=FALSE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           xlab = "Relative Feature importance",lwd.ci = 1.4,boxsize = 0.1)
