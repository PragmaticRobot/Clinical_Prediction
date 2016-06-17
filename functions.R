####################################
######### 99 - Functions ###########
####################################

cv_RF <- function(feat,out)
{
  RF_model <- foreach(ntree = rep(10000,4), .combine = combine, .multicombine=TRUE, 
                      .packages = "randomForest") %dopar% randomForest(feat,out,ntree=ntree)
  imp <- sort(RF_model$importance, partial=NULL,decreasing=TRUE, index.return=TRUE)
  predict <- predict(RF_model, feat, type = "response" , predict.all=FALSE)
  GRE <- list(predict=predict, Index=imp$ix, imps=imp$x)
  return(GRE)
}

cv_mod <- function(feat, out, stand, isWO)
{
  if(stand == 1)
  {
    lassoMOD <- cv.glmnet(x=as.matrix(feat),y=out, standardize = TRUE)
  } else {
    lassoMOD <- cv.glmnet(x=as.matrix(feat),y=out, standardize = FALSE)
  }
  cutoff <- min(lassoMOD$cvup)
  means <- lassoMOD$cvm
  if(isWO == 1)
  {
    ind <- min(which(means<cutoff)) # now we know which model won
  } else {
    ind <- max(which(means<cutoff)) # now we know which model won
  }
  predict <- predict(lassoMOD, newx=as.matrix(feat), s= lassoMOD$lambda[ind])
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
  x1 <- ceiling(max(sorted_list, na.rm = TRUE))
  x2 <- floor(min(sorted_list, na.rm = TRUE))
  yy <- rep(N, length(gg))+ offs
  plot(gg,yy,col='blue',pch=16,ylim = c(1, N),
       xlim = c(x2, x1),cex=0.3,yaxt='n', main=tit,
       xlab="",ylab="",cex.lab=0.7)
  # mean segment
  xx1 <- fivenum(gg,na.rm = TRUE)
  xq1 <- xx1[2];  xq2 <- xx1[3]; xq3 <- xx1[4]
  segments(xq2-0.05, N-0.5, xq2-0.05, N+0.5, col = 'red',lwd = 2) # center (median)
  segments(xq1-0.05, N-0.5, xq1-0.05, N+0.5, col = 'forestgreen',lwd = 2) # left (Q1)
  segments(xq3-0.05, N-0.5, xq3-0.05, N+0.5, col = 'forestgreen',lwd = 2) # right (Q3)
  segments(xq1-0.05, N+0.5, xq3-0.05, N+0.5, col = 'forestgreen',lwd = 2) # top
  segments(xq1-0.05, N-0.5, xq3-0.05, N-0.5, col = 'forestgreen',lwd = 2) # bottom
  axis(side = 2, at = N,paste(nam))
  for (i in 2:N){
    x <- sorted_list[i,]
    nam <- rownames(sorted_list)[i]
    x <- x[!is.na(x)]
    offs <- runif(length(x),-0.35,0.35)
    gg <- x
    yy <- rep(N - (i-1), length(gg))+ offs
    points(gg,yy,col='blue',pch=16,cex=0.3)
    yc <- N - (i-1) # the center of the y-axis for this feature
    xx1 <- fivenum(gg,na.rm = TRUE)
    xq1 <- xx1[2];  xq2 <- xx1[3]; xq3 <- xx1[4]
    segments(xq2-0.05, yc-0.5, xq2-0.05, yc+0.5, col = 'red',lwd = 2) # center (median)
    segments(xq1-0.05, yc-0.5, xq1-0.05, yc+0.5, col = 'forestgreen',lwd = 2) # left (Q1)
    segments(xq3-0.05, yc-0.5, xq3-0.05, yc+0.5, col = 'forestgreen',lwd = 2) # right (Q3)
    segments(xq1-0.05, yc+0.5, xq3-0.05, yc+0.5, col = 'forestgreen',lwd = 2) # top
    segments(xq1-0.05, yc-0.5, xq3-0.05, yc-0.5, col = 'forestgreen',lwd = 2) # bottom
    axis(side = 2, at = yc,paste(nam))
  }
}
