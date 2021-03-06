####################################
######### 99 - Functions ###########
####################################
SecondOrderThis <- function(FirstOrder)
{
  SecondOrder <- FirstOrder
  names2 <- c(colnames(FirstOrder))
  for (i in c(1:dim(FirstOrder)[2])) {
    for (j in c(i:dim(FirstOrder)[2])){
      if (is.factor(FirstOrder[,i]) & is.factor(FirstOrder[,j]) & identical(FirstOrder[,i],FirstOrder[,j])) # If they're the same and both factors, skip
      {
        next
      } else if (is.factor(FirstOrder[,i]) & is.factor(FirstOrder[,j]) & !identical(FirstOrder[,i],FirstOrder[,j])) # If not the same and both factors, interact
      {
        name_new <- paste(colnames(FirstOrder)[i], ':',colnames(FirstOrder)[j], sep="")
        trash <- interaction(FirstOrder[,i],FirstOrder[,j],sep = ":")
        df_new <- data.frame(new_name = (trash))
        colnames(df_new) <- name_new
        # now check for df3, whether any of the "cells" in the table are empty
        if (min(table(df_new))==0){
          next
        } else {
          SecondOrder <- cbind(SecondOrder,(df_new))
          names2 <- c(names2,name_new)
        }
      } else if ((is.factor(FirstOrder[,i]) & !is.factor(FirstOrder[,j])) | 
                 (!is.factor(FirstOrder[,i]) & is.factor(FirstOrder[,j])))
      {
        name_new <- paste(colnames(FirstOrder)[i], ':',colnames(FirstOrder)[j], sep="")
        trash <- interaction(FirstOrder[,i],FirstOrder[,j],sep = ":")
        df_new <- data.frame(new_name = (trash))
        colnames(df_new) <- name_new
        # now check for df3, whether any of the "cells" in the table are empty
        SecondOrder <- cbind(SecondOrder,(df_new))
        names2 <- c(names2,name_new)
      } else
      {
        name_new <- paste(colnames(FirstOrder)[i], ':',colnames(FirstOrder)[j], sep="")
        trash <- (FirstOrder[,i]*FirstOrder[,j])
        df_new <- data.frame(name_new = (trash))
        colnames(df_new) <- name_new
        # now check for df3, whether any of the "cells" in the table are empty
        SecondOrder <- cbind(SecondOrder,(df_new))
        names2 <- c(names2,name_new)
      }
    }
  }
  colnames(SecondOrder) <- c(names2)
  return(SecondOrder)
}

cv_RF <- function(feat,out)
{
  RF_model <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% randomForest(feat,out,ntree=ntree, corr.bias=FALSE)
  imp <- sort(RF_model$importance, partial=NULL,decreasing=TRUE, index.return=TRUE)
  predict <- predict(RF_model, feat, type = "response" , predict.all=FALSE)
  GRE <- list(predict=predict, Index=imp$ix, imps=imp$x)
  return(GRE)
}

rf_crossvd <- function(feat,out)
{
  RF_model <- foreach(ntree = rep(12500,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% rfcv(feat,out,ntree=ntree, cv.fold=4)
  # imp <- sort(RF_model$importance, partial=NULL,decreasing=TRUE, index.return=TRUE)
  # predict <- predict(RF_model, feat, type = "response" , predict.all=FALSE)
  # GRE <- list(predict=predict, Index=imp$ix, imps=imp$x)
  return(RF_model)
}

cv_mod <- function(feat, out, stand, isWO, nf)
{
  if(stand == 1)
  {
    lassoMOD <- cv.glmnet(data.matrix(feat),out, standardize = TRUE, alpha = 1, nlambda = 1000, nfolds = nf)
  } else {
    lassoMOD <- cv.glmnet(data.matrix(feat),out, standardize = FALSE, alpha = 1, nlambda = 1000, nfolds = nf)
  }
  cutoff <- min(lassoMOD$cvup)
  means <- lassoMOD$cvm
  if(isWO == 1)
  {
    ind <- max(which(means<cutoff)) # now we know which model won
  } else {
    ind <- max(which(means<cutoff)) # now we know which model won
  }
  predict <- predict(lassoMOD, newx=data.matrix(feat), s= lassoMOD$lambda[ind])
  cc <- coef(lassoMOD$glmnet.fit, s = lassoMOD$lambda[ind])
  summ <- summary(cc)
  GRE <- list(predict=predict, Index = summ$i, coefs = summ$x)
  return(GRE)
}

## for the cleanup functions flag:
## 1: predictions
## 2: indices of winning variables (LASSO) or order of variable importance (RF)
## 3: value for coefficients of winning variables, or importance of variable (RF)
cleanup_LASSO <- function(pred1)
{
  predictions1 <- pred1[[1]]$predict
  inds1 <- pred1[[1]]$Index
  coefs1 <- pred1[[1]]$coefs
  for (i in 2:100)
  {
    predictions1 <- cbind(predictions1,pred1[[i]]$predict)
    inds1 <- qpcR:::cbind.na(inds1,pred1[[i]]$Index)
    coefs1 <- qpcR:::cbind.na(coefs1,pred1[[i]]$coefs)
  }
  LASSO_pred1 <- list(Predictions=predictions1,Indices=inds1,Coeffs=coefs1)
  return(LASSO_pred1)
}

cleanup_RF <- function(pred1)
{
  predictions1 <- pred1[[1]]$predict
  inds1 <- pred1[[1]]$Index
  imps1 <- pred1[[1]]$imps
  for (i in 2:100)
  {
    predictions1 <- cbind(predictions1,pred1[[i]]$predict)
    inds1 <- qpcR:::cbind.na(inds1,pred1[[i]]$Index)
    imps1 <- qpcR:::cbind.na(imps1,pred1[[i]]$imps)
  }
  RF_pred1 <- list(Predictions=predictions1,Indices=inds1,Imps=imps1)
  return(RF_pred1)
}

# plot_sideways plots the variable importance for the cross-validation of
# the random forest algorithm, it takes the "sorted list", a data frame with
# representing features and columns representing cross-validations, the rows
# should have the feature name as rownames
#
# The second input is N, if N is zero, the function uses the number of features
# in the list as N, otherwise it uses the value provided when the function is
# called
#
# The third input is a title string, used for the "main" title for the plot
#
#
plot_sideways <- function(sorted_list,N,tit)
{
  par(las=1)
  par(mar=c(4,20,4,2))
  x <- sorted_list[1,]
  nam <- rownames(sorted_list)[1]
  x <- x[!is.na(x)]
  offs <- runif(length(x),-0.35,0.35)
  gg <- x
  if(N == 0){
    N <- dim(sorted_list)[1]
  }
  xg1 <- max(sorted_list, na.rm = TRUE)
  xg2 <- min(sorted_list, na.rm = TRUE)
  if (xg1>1) {
    x1 <- ceiling(xg1)
  }
  else {
    x1 <- round(xg1,2)
  }
  if (xg2>1){
    x2 <- floor(xg2)
  }
  else {
    x2 <- round(xg2,2)
  }
  yy <- rep(N, length(gg))+ offs
  plot(gg,yy,col='darkcyan',pch=16,ylim = c(1, N),
       xlim = c(x2, x1),cex=0.3,yaxt='n', main=tit,
       xlab="",ylab="",cex.lab=0.7)
  # mean segment
  xx1 <- fivenum(gg,na.rm = TRUE)
  xq1 <- xx1[2];  xq2 <- xx1[3]; xq3 <- xx1[4]
  segments(xq2, N-0.5, xq2, N+0.5, col = 'firebrick2',lwd = 2) # center (median)
  segments(xq1, N-0.5, xq1, N+0.5, col = 'forestgreen',lwd = 2) # left (Q1)
  segments(xq3, N-0.5, xq3, N+0.5, col = 'forestgreen',lwd = 2) # right (Q3)
  segments(xq1, N+0.5, xq3, N+0.5, col = 'forestgreen',lwd = 2) # top
  segments(xq1, N-0.5, xq3, N-0.5, col = 'forestgreen',lwd = 2) # bottom
  axis(side = 2, at = N,paste(nam))
  for (i in 2:N){
    x <- sorted_list[i,]
    nam <- rownames(sorted_list)[i]
    x <- x[!is.na(x)]
    offs <- runif(length(x),-0.35,0.35)
    gg <- x
    yy <- rep(N - (i-1), length(gg))+ offs
    points(gg,yy,col='darkcyan',pch=16,cex=0.3)
    yc <- N - (i-1) # the center of the y-axis for this feature
    xx1 <- fivenum(gg,na.rm = TRUE)
    xq1 <- xx1[2];  xq2 <- xx1[3]; xq3 <- xx1[4]
    segments(xq2, yc-0.5, xq2, yc+0.5, col = 'firebrick2',lwd = 2) # center (median)
    segments(xq1, yc-0.5, xq1, yc+0.5, col = 'forestgreen',lwd = 2) # left (Q1)
    segments(xq3, yc-0.5, xq3, yc+0.5, col = 'forestgreen',lwd = 2) # right (Q3)
    segments(xq1, yc+0.5, xq3, yc+0.5, col = 'forestgreen',lwd = 2) # top
    segments(xq1, yc-0.5, xq3, yc-0.5, col = 'forestgreen',lwd = 2) # bottom
    axis(side = 2, at = yc,paste(nam))
  }
  segments(x0=0,y0=0,x1=0,y1=N+1, col = "black",lwd=1.5, lty=2)
}

# The plotting function is specific here, it has to take either 2 inputs or 4 inputs
# for 2 inputs: Linear LASSO and Linear RF
# for 4 inputs: Linear LASSO, quad LASSO, Linear RF, quad RF
# ORDER IS IMPORTANT!!!
# modlist should have the models listed in exactly the same order
plot_LSRF_fits <- function(allins,modlist,relative,out,outname,measure)
{
  if (ncol(allins)==2){
    if (relative == 1){
      plot(allins[,1],out, xlim = c(round(min(out),1), round(max(out),1)), 
           ylim = c(round(min(out),1), round(max(out),1)), col='firebrick2',
           ylab=paste('Actual Relative Change in',outname),
           xlab=paste('Predicted Relative Change in',outname),
           cex.lab=0.8,pch=21, cex=0.8,lwd=2)
    } else {
      plot(allins[,1],out, xlim = c(round(min(out),1), round(max(out),1)), 
           ylim = c(round(min(out),1), round(max(out),1)), col='firebrick2',
           ylab=paste('Actual Change in',outname),
           xlab=paste('Predicted Change in',outname),
           cex.lab=0.8,pch=21, cex=0.8,lwd=2)
    }
    abline(modlist[[1]], col='firebrick2',lwd=2)
    points(allins[,2],out, col = 'forestgreen',pch=22, cex=0.8,lwd=2)
    abline(modlist[[2]], col='forestgreen',lwd=2)
    abline(0,1, col='black',pch=16,cex=0.6,lwd=2)
    if (measure == 1){
      r1 <- round(sqrt(mean((out-allins[,1])^2)), digits = 2)
      r2 <- round(sqrt(mean((out-allins[,2])^2)), digits = 2)
      meString <- 'RMSE'
    } else if (measure == 2){
      blah <- summary(modlist[[1]],correlation=TRUE)
      r1 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[2]],correlation=TRUE)
      r2 <- round(blah$correlation, digits = 2)
      meString <- 'Corr'
    } else if (measure == 3){
      blah <- summary(modlist[[1]])
      r1 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[2]])
      r2 <- round(blah$adj.r.squared, digits = 2)
      meString <- 'R^2'
    } else {
      message('unknown outcome measure')
    }
    legend("topleft",c(paste('Linear LASSO ',meString ,' = ',as.character(r1)),
                       paste('Linear Random Forests ',meString,' = ',as.character(r2))),
           cex=0.8, col=c('firebrick2','forestgreen'),
           lwd=c(2.5,2.5),lty=c(1,1),bty="n")
  } else if (ncol(allins)==4) {
    if (relative==1){
      plot(allins[,1],out, xlim = c(round(min(out),1), round(max(out),1)), 
           ylim = c(round(min(out),1), round(max(out),1)), col='firebrick2',
           ylab=paste('Actual Relative Change in',outname),
           xlab=paste('Predicted Relative Change in',outname),
           cex.lab=0.8,pch=21, cex=0.8,lwd=2)
    } else {
      plot(allins[,1],out, xlim = c(round(min(out),1), round(max(out),1)), 
           ylim = c(round(min(out),1), round(max(out),1)), col='firebrick2',
           ylab=paste('Actual Change in',outname),
           xlab=paste('Predicted Change in',outname),
           cex.lab=0.8,pch=21, cex=0.8,lwd=2)
    }
    abline(modlist[[1]], col='firebrick2')
    points(allins[,2],out, col = 'darkcyan',pch=22, cex=0.8)
    abline(modlist[[2]], col='darkcyan')
    points(allins[,3],out, col = 'forestgreen',pch=24, cex=0.8)
    abline(modlist[[3]], col='forestgreen')
    points(allins[,4],out, col = 'darkorange1',pch=25, cex=0.8)
    abline(modlist[[4]], col='darkorange1')
    abline(0,1, col='black',pch=16,cex=0.6,lwd=2)
    if (measure == 1){
      r1 <- round(sqrt(mean((out-allins[,1])^2)), digits = 2)
      r2 <- round(sqrt(mean((out-allins[,2])^2)), digits = 2)
      r3 <- round(sqrt(mean((out-allins[,3])^2)), digits = 2)
      r4 <- round(sqrt(mean((out-allins[,4])^2)), digits = 2)
      meString <- 'RMSE'
    } else if (measure == 2){
      blah <- summary(modlist[[1]],correlation=TRUE)
      r1 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[2]],correlation=TRUE)
      r2 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[3]],correlation=TRUE)
      r3 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[4]],correlation=TRUE)
      r4 <- round(blah$correlation, digits = 2)
      meString <- 'Corr'
    } else if (measure == 3){
      blah <- summary(modlist[[1]])
      r1 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[2]])
      r2 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[3]])
      r3 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[4]])
      r4 <- round(blah$adj.r.squared, digits = 2)
      meString <- 'R^2'
    } else {
      message('unknown outcome measure')
    }
    legend("topleft",c(paste('Linear LASSO ',meString,' = ',as.character(r1)),
                       paste('Quadratic LASSO ',meString,' = ',as.character(r2)),
                       paste('Linear Random Forests ',meString,' = ',as.character(r3)),
                       paste('Quadratic Random Forests ',meString,' = ',as.character(r4))),
           cex=0.8, col=c('firebrick2','darkcyan','forestgreen','darkorange1'),
           lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5),bty="n")
  } else if (ncol(allins)==6) {
    if (relative==1){
      plot(allins[,1],out, xlim = c(round(min(out),1), round(max(out),1)), 
           ylim = c(round(min(out),1), round(max(out),1)), col='firebrick2',
           ylab=paste('Actual Relative Change in',outname),
           xlab=paste('Predicted Relative Change in',outname),
           cex.lab=0.8,pch=21, cex=0.8,lwd=2)
    } else {
      plot(allins[,1],out, xlim = c(round(min(out),1), round(max(out),1)), 
           ylim = c(round(min(out),1), round(max(out),1)), col='firebrick2',
           ylab=paste('Actual Change in',outname),
           xlab=paste('Predicted Change in',outname),
           cex.lab=0.8,pch=21, cex=0.8,lwd=2)
    }
    abline(modlist[[1]], col='firebrick2',lwd=3)
    points(allins[,2],out, col = 'darkcyan',pch=22, cex=0.8)
    abline(modlist[[2]], col='darkcyan',lwd=3)
    points(allins[,3],out, col = 'forestgreen',pch=24, cex=0.8)
    abline(modlist[[3]], col='forestgreen',lwd=3)
    points(allins[,4],out, col = 'darkorange1',pch=25, cex=0.8)
    abline(modlist[[4]], col='darkorange1',lwd=3)
    points(allins[,5],out, col = 'khaki4',pch=3, cex=0.8)
    abline(modlist[[5]], col='khaki4',lwd=3)
    points(allins[,6],out, col = 'slateblue2',pch=8, cex=0.8)
    abline(modlist[[6]], col='slateblue2',lwd=3)
    abline(0,1, col='black',pch=16,cex=0.6,lwd=2)
    if (measure == 1){
      r1 <- round(sqrt(mean((out-allins[,1])^2)), digits = 2)
      r2 <- round(sqrt(mean((out-allins[,2])^2)), digits = 2)
      r3 <- round(sqrt(mean((out-allins[,3])^2)), digits = 2)
      r4 <- round(sqrt(mean((out-allins[,4])^2)), digits = 2)
      r5 <- round(sqrt(mean((out-allins[,5])^2)), digits = 2)
      r6 <- round(sqrt(mean((out-allins[,6])^2)), digits = 2)
      meString <- 'RMSE'
    } else if (measure == 2){
      blah <- summary(modlist[[1]],correlation=TRUE)
      r1 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[2]],correlation=TRUE)
      r2 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[3]],correlation=TRUE)
      r3 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[4]],correlation=TRUE)
      r4 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[5]],correlation=TRUE)
      r5 <- round(blah$correlation, digits = 2)
      blah <- summary(modlist[[6]],correlation=TRUE)
      r6 <- round(blah$correlation, digits = 2)
      meString <- 'Corr'
    } else if (measure == 3){
      blah <- summary(modlist[[1]])
      r1 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[2]])
      r2 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[3]])
      r3 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[4]])
      r4 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[5]])
      r5 <- round(blah$adj.r.squared, digits = 2)
      blah <- summary(modlist[[6]])
      r6 <- round(blah$adj.r.squared, digits = 2)
      meString <- 'R^2'
    } else {
      message('unknown outcome measure')
    }
    legend("topleft",c(paste('Linear LASSO ',meString,' = ',as.character(r1)),
                       paste('Quadratic LASSO ',meString,' = ',as.character(r2)),
                       paste('Linear Random Forests ',meString,' = ',as.character(r3)),
                       paste('Quadratic Random Forests ',meString,' = ',as.character(r4)),
                       paste('NoExtrap LASSO ',meString,' = ',as.character(r5)),
                       paste('NoExtrap Random Forests ',meString,' = ',as.character(r6))),
           cex=0.8, col=c('firebrick2','darkcyan','forestgreen','darkorange1','khaki4','slateblue2'),
           lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),bty="n")
  } else {
    message("Unsupported input size")
  }
}

## Add a function DoHists here, have it return skewness
## Bar plots showing quality of model fits
# Columns: 1- Linear LASSO, 2- Quad LASSO, 3- Linear RF, 4- Quad RF
DoHists <- function(allins,out,outname)
{
  if (ncol(allins)==2){
    difFM <- matrix(, nrow = 26, ncol = 2)
    difFM[,1] <- abs(allins[,1]-out)
    difFM[,2] <- abs(allins[,2]-out)
    bks <- ceiling(apply(difFM,2,max,na.rm = TRUE)) # find the xlim of each hist
    D1 <- hist(difFM[,1], breaks = bks[1], xlim = c(0,bks[1]+1), ylim = c(0,12),plot=FALSE)
    D2 <- hist(difFM[,2], breaks = bks[2], xlim = c(0,bks[2]+1), ylim = c(0,12),plot=FALSE)
    # Grouped Bar Plot
    counts <- matrix(, nrow = max(bks), ncol = 2)
    counts[,1] <- padzeros(D1$counts,max(bks)-bks[1],"right")
    counts[,2] <- padzeros(D2$counts,max(bks)-bks[2],"right")
    histData <- matrix(counts,ncol = max(bks), byrow=T)
    rownames(histData) <- c("Linear LASSO","Linear Random Forests")
    colns <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12")
    colnames(histData) <- colns[1:max(bks)]
    barplot(histData, main=paste("Prediction Accuracy for",outname),
            xlab=paste("Absolute prediction error (in units of)",outname,")"),
            ylab = "number of subjects",
            col=c("firebrick2","forestgreen"),
            legend=rownames(histData),beside=TRUE,ylim = c(0,max(counts)+1))
    sks <- apply(counts,2,e1071::skewness,na.rm=TRUE)
  } else if (ncol(allins)==4){
    difFM <- matrix(, nrow = 26, ncol = 4)
    difFM[,1] <- abs(allins[,1]-out)
    difFM[,2] <- abs(allins[,2]-out)
    difFM[,3] <- abs(allins[,3]-out)
    difFM[,4] <- abs(allins[,4]-out)
    bks <- ceiling(apply(difFM,2,max,na.rm = TRUE)) # find the xlim of each hist
    D1 <- hist(difFM[,1], breaks = c(0:bks[1]), xlim = c(0,bks[1]+1), ylim = c(0,12),plot=FALSE)
    D2 <- hist(difFM[,2], breaks = c(0:bks[2]), xlim = c(0,bks[2]+1), ylim = c(0,12),plot=FALSE)
    D3 <- hist(difFM[,3], breaks = c(0:bks[3]), xlim = c(0,bks[3]+1), ylim = c(0,12),plot=FALSE)
    D4 <- hist(difFM[,4], breaks = c(0:bks[4]), xlim = c(0,bks[4]+1), ylim = c(0,12),plot=FALSE)
    # Grouped Bar Plot
    counts <- matrix(, nrow = max(bks), ncol = 4)
    counts[,1] <- padzeros(D1$counts,max(bks)-bks[1],"right")
    counts[,2] <- padzeros(D2$counts,max(bks)-bks[2],"right")
    counts[,3] <- padzeros(D3$counts,max(bks)-bks[3],"right")
    counts[,4] <- padzeros(D4$counts,max(bks)-bks[4],"right")
    histData <- matrix(counts,ncol = max(bks), byrow=T)
    rownames(histData) <- c("Linear LASSO","Quadratic LASSO","Linear Random Forests","Quadratic RandomForests")
    colns <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12")
    colnames(histData) <- colns[1:max(bks)]
    barplot(histData, main=paste("Prediction Accuracy for",outname),
            xlab=paste("Absolute prediction error (in units of)",outname,")"),
            ylab = "number of subjects",
            col=c("firebrick2","darkcyan","forestgreen","darkorange1"),
            legend=rownames(histData),beside=TRUE,ylim = c(0,max(counts)+1))
    sks <- apply(counts,2,e1071::skewness,na.rm=TRUE)
  } else if (ncol(allins)==6){
    difFM <- matrix(, nrow = 26, ncol = 6)
    difFM[,1] <- abs(allins[,1]-out)
    difFM[,2] <- abs(allins[,2]-out)
    difFM[,3] <- abs(allins[,3]-out)
    difFM[,4] <- abs(allins[,4]-out)
    difFM[,5] <- abs(allins[,5]-out)
    difFM[,6] <- abs(allins[,6]-out)
    bks <- ceiling(apply(difFM,2,max,na.rm = TRUE)) # find the xlim of each hist
    D1 <- hist(difFM[,1], breaks = c(0:bks[1]), xlim = c(0,bks[1]+1), ylim = c(0,12),plot=FALSE)
    D2 <- hist(difFM[,2], breaks = c(0:bks[2]), xlim = c(0,bks[2]+1), ylim = c(0,12),plot=FALSE)
    D3 <- hist(difFM[,3], breaks = c(0:bks[3]), xlim = c(0,bks[3]+1), ylim = c(0,12),plot=FALSE)
    D4 <- hist(difFM[,4], breaks = c(0:bks[4]), xlim = c(0,bks[4]+1), ylim = c(0,12),plot=FALSE)
    D5 <- hist(difFM[,5], breaks = c(0:bks[5]), xlim = c(0,bks[5]+1), ylim = c(0,12),plot=FALSE)
    D6 <- hist(difFM[,6], breaks = c(0:bks[6]), xlim = c(0,bks[6]+1), ylim = c(0,12),plot=FALSE)
    # Grouped Bar Plot
    counts <- matrix(, nrow = max(bks), ncol = 6)
    counts[,1] <- padzeros(D1$counts,max(bks)-bks[1],"right")
    counts[,2] <- padzeros(D2$counts,max(bks)-bks[2],"right")
    counts[,3] <- padzeros(D3$counts,max(bks)-bks[3],"right")
    counts[,4] <- padzeros(D4$counts,max(bks)-bks[4],"right")
    counts[,5] <- padzeros(D5$counts,max(bks)-bks[5],"right")
    counts[,6] <- padzeros(D6$counts,max(bks)-bks[6],"right")
    histData <- matrix(counts,ncol = max(bks), byrow=T)
    rownames(histData) <- c("Linear LASSO","Quadratic LASSO",
                            "Linear Random Forests","Quadratic RandomForests",
                            "NoExtrap LASSO","NoExtrap Random Forests")
    colns <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12")
    colnames(histData) <- colns[1:max(bks)]
    barplot(histData, main=paste("Prediction Accuracy for",outname),
            xlab=paste("Absolute prediction error (in units of)",outname,")"),
            ylab = "number of subjects",
            col=c("firebrick2","darkcyan","forestgreen","darkorange1","khaki4","slateblue2"),
            legend=rownames(histData),beside=TRUE,ylim = c(0,max(counts)+1))
    sks <- apply(counts,2,e1071::skewness,na.rm=TRUE)
  } else {
    message("Unsupported input size")
  }
  return(sks)
}

# Here, rewrite DoHists to make line plot, everything is the same
# except for the barplot command
DoLines <- function(allins,out,outname)
{
  if (ncol(allins)==2){
    difFM <- matrix(, nrow = 26, ncol = 2)
    difFM[,1] <- abs(allins[,1]-out)
    difFM[,2] <- abs(allins[,2]-out)
    bks <- ceiling(apply(difFM,2,max,na.rm = TRUE)) # find the xlim of each hist
    D1 <- hist(difFM[,1], breaks = bks[1], xlim = c(0,bks[1]+1), ylim = c(0,12),plot=FALSE)
    D2 <- hist(difFM[,2], breaks = bks[2], xlim = c(0,bks[2]+1), ylim = c(0,12),plot=FALSE)
    # Grouped Bar Plot
    counts <- matrix(, nrow = max(bks), ncol = 2)
    counts[,1] <- padzeros(D1$counts,max(bks)-bks[1],"right")
    counts[,2] <- padzeros(D2$counts,max(bks)-bks[2],"right")
    histData <- matrix(counts,ncol = max(bks), byrow=T)
    rownames(histData) <- c("Linear LASSO","Linear Random Forests")
    colns <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12")
    colnames(histData) <- colns[1:max(bks)]
    xx1 <- seq(0.5,ncol(histData)-0.5,by = 1)
    plot(xx1,histData[1,], pch=21, lwd=2.5,col="firebrick2", 
         xlim=c(0,ncol(histData)),ylim=c(0,max(histData,na.rm = TRUE)),
         xaxt="n",xlab=paste("Absolute prediction error (in units of",outname,")"),
         main=paste("Prediction Accuracy for",outname),
         ylab = "number of subjects")
    lines(xx1,histData[1,],lwd=2.5,lty=1,col="firebrick2")
    points(xx1,histData[2,],pch=24,lwd=2.5,col="forestgreen")
    lines(xx1,histData[2,],lwd=2.5,lty=1,col="forestgreen")
    axis(side=1,at=xx1,labels=colnames(histData))
    sks <- apply(counts,2,e1071::skewness,na.rm=TRUE)
    sks2 <- round(sks,digits = 2)
    legend("topright",c(paste(rownames(histData)[1],", skew = ",as.character(sks2[1])),
                        paste(rownames(histData)[2],", skew = ",as.character(sks2[2]))),
           cex = 1, col=c('firebrick2','forestgreen'),
           lty=c(1,1),lwd=c(2.5,2.5),bty="n")
  } else if (ncol(allins)==4){
    difFM <- matrix(, nrow = 26, ncol = 4)
    difFM[,1] <- abs(allins[,1]-out)
    difFM[,2] <- abs(allins[,2]-out)
    difFM[,3] <- abs(allins[,3]-out)
    difFM[,4] <- abs(allins[,4]-out)
    bks <- ceiling(apply(difFM,2,max,na.rm = TRUE)) # find the xlim of each hist
    D1 <- hist(difFM[,1], breaks = c(0:bks[1]), xlim = c(0,bks[1]+1), ylim = c(0,12),plot=FALSE)
    D2 <- hist(difFM[,2], breaks = c(0:bks[2]), xlim = c(0,bks[2]+1), ylim = c(0,12),plot=FALSE)
    D3 <- hist(difFM[,3], breaks = c(0:bks[3]), xlim = c(0,bks[3]+1), ylim = c(0,12),plot=FALSE)
    D4 <- hist(difFM[,4], breaks = c(0:bks[4]), xlim = c(0,bks[4]+1), ylim = c(0,12),plot=FALSE)
    # Line Plot
    counts <- matrix(, nrow = max(bks), ncol = 4)
    counts[,1] <- padzeros(D1$counts,max(bks)-bks[1],"right")
    counts[,2] <- padzeros(D2$counts,max(bks)-bks[2],"right")
    counts[,3] <- padzeros(D3$counts,max(bks)-bks[3],"right")
    counts[,4] <- padzeros(D4$counts,max(bks)-bks[4],"right")
    histData <- matrix(counts,ncol = max(bks), byrow=T)
    rownames(histData) <- c("Linear LASSO","Quadratic LASSO","Linear Random Forests","Quadratic RandomForests")
    colns <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12")
    colnames(histData) <- colns[1:max(bks)]
    xx1 <- seq(0.5,ncol(histData)-0.5,by = 1)
    plot(xx1,histData[1,], pch=21, lwd=2.5,col="firebrick2", 
         xlim=c(0,ncol(histData)),ylim=c(0,max(histData,na.rm = TRUE)),
         xaxt="n",xlab=paste("Absolute prediction error (in units of",outname,")"),
         main=paste("Prediction Accuracy for",outname),
         ylab = "number of subjects")
    lines(xx1,histData[1,],lwd=2.5,lty=1,col="firebrick2")
    points(xx1,histData[2,],pch=22,lwd=2.5,col="darkcyan")
    lines(xx1,histData[2,],lwd=2.5,lty=1,col="darkcyan")
    points(xx1,histData[3,],pch=24,lwd=2.5,col="forestgreen")
    lines(xx1,histData[3,],lwd=2.5,lty=1,col="forestgreen")
    points(xx1,histData[4,],pch=25,lwd=2.5,col="darkorange1")
    lines(xx1,histData[4,],lwd=2.5,lty=1,col="darkorange1")
    axis(side=1,at=xx1,labels=colnames(histData))
    sks <- apply(counts,2,e1071::skewness,na.rm=TRUE)
    sks2 <- round(sks,digits = 2)
    legend("topright",c(paste(rownames(histData)[1],", skew = ",as.character(sks2[1])),
                        paste(rownames(histData)[2],", skew = ",as.character(sks2[2])),
                        paste(rownames(histData)[3],", skew = ",as.character(sks2[3])),
                        paste(rownames(histData)[4],", skew = ",as.character(sks2[4]))),
           cex = 1, col=c('firebrick2','darkcyan','forestgreen','darkorange1'),
           lty=c(1,1,1,1),lwd=c(3,3,3,3),bty="n")
  } else if (ncol(allins)==6){
    difFM <- matrix(, nrow = 26, ncol = 6)
    difFM[,1] <- abs(allins[,1]-out)
    difFM[,2] <- abs(allins[,2]-out)
    difFM[,3] <- abs(allins[,3]-out)
    difFM[,4] <- abs(allins[,4]-out)
    difFM[,5] <- abs(allins[,5]-out)
    difFM[,6] <- abs(allins[,6]-out)
    bks <- ceiling(apply(difFM,2,max,na.rm = TRUE)) # find the xlim of each hist
    D1 <- hist(difFM[,1], breaks = c(0:bks[1]), xlim = c(0,bks[1]+1), ylim = c(0,12),plot=FALSE)
    D2 <- hist(difFM[,2], breaks = c(0:bks[2]), xlim = c(0,bks[2]+1), ylim = c(0,12),plot=FALSE)
    D3 <- hist(difFM[,3], breaks = c(0:bks[3]), xlim = c(0,bks[3]+1), ylim = c(0,12),plot=FALSE)
    D4 <- hist(difFM[,4], breaks = c(0:bks[4]), xlim = c(0,bks[4]+1), ylim = c(0,12),plot=FALSE)
    D5 <- hist(difFM[,5], breaks = c(0:bks[5]), xlim = c(0,bks[5]+1), ylim = c(0,12),plot=FALSE)
    D6 <- hist(difFM[,6], breaks = c(0:bks[6]), xlim = c(0,bks[6]+1), ylim = c(0,12),plot=FALSE)
    # Line Plot
    counts <- matrix(, nrow = max(bks), ncol = 6)
    counts[,1] <- padzeros(D1$counts,max(bks)-bks[1],"right")
    counts[,2] <- padzeros(D2$counts,max(bks)-bks[2],"right")
    counts[,3] <- padzeros(D3$counts,max(bks)-bks[3],"right")
    counts[,4] <- padzeros(D4$counts,max(bks)-bks[4],"right")
    counts[,5] <- padzeros(D5$counts,max(bks)-bks[5],"right")
    counts[,6] <- padzeros(D6$counts,max(bks)-bks[6],"right")
    histData <- matrix(counts,ncol = max(bks), byrow=T)
    rownames(histData) <- c("Linear LASSO","Quadratic LASSO",
                            "Linear Random Forests","Quadratic RandomForests",
                            "NoExtrap LASSO","NoExtrap Random Forests")
    colns <- c("0-1","1-2","2-3","3-4","4-5","5-6","6-7","7-8","8-9","9-10","10-11","11-12")
    colnames(histData) <- colns[1:max(bks)]
    xx1 <- seq(0.5,ncol(histData)-0.5,by = 1)
    plot(xx1,histData[1,], pch=21, lwd=3,col="firebrick2", 
         xlim=c(0,ncol(histData)),ylim=c(0,max(histData,na.rm = TRUE)),
         xaxt="n",xlab=paste("Absolute prediction error (in units of",outname,")"),
         main=paste("Prediction Accuracy for",outname),
         ylab = "number of subjects")
    lines(xx1,histData[1,],lwd=3,lty=1,col="firebrick2")
    points(xx1,histData[2,],pch=22,lwd=3,col="darkcyan")
    lines(xx1,histData[2,],lwd=3,lty=1,col="darkcyan")
    points(xx1,histData[3,],pch=24,lwd=3,col="forestgreen")
    lines(xx1,histData[3,],lwd=3,lty=1,col="forestgreen")
    points(xx1,histData[4,],pch=25,lwd=3,col="darkorange1")
    lines(xx1,histData[4,],lwd=3,lty=1,col="darkorange1")
    points(xx1,histData[5,],pch=3,lwd=3,col="khaki4")
    lines(xx1,histData[5,],lwd=3,lty=1,col="khaki4")
    points(xx1,histData[6,],pch=8,lwd=3,col="slateblue2")
    lines(xx1,histData[6,],lwd=3,lty=1,col="slateblue2")
    axis(side=1,at=xx1,labels=colnames(histData))
    sks <- apply(counts,2,e1071::skewness,na.rm=TRUE)
    sks2 <- round(sks,digits = 2)
    legend("topright",c(paste(rownames(histData)[1],", skew = ",as.character(sks2[1])),
                        paste(rownames(histData)[2],", skew = ",as.character(sks2[2])),
                        paste(rownames(histData)[3],", skew = ",as.character(sks2[3])),
                        paste(rownames(histData)[4],", skew = ",as.character(sks2[4])),
                        paste(rownames(histData)[5],", skew = ",as.character(sks2[5])),
                        paste(rownames(histData)[6],", skew = ",as.character(sks2[6]))),
           cex = 1, col=c('firebrick2','darkcyan','forestgreen','darkorange1','khaki4','slateblue2'),
           lty=c(1,1,1,1,1,1),lwd=c(3,3,3,3,3,3),bty="n")
  } else {
    message("Unsupported input size")
  }
  return(sks)
}