
## This script will plot the results of the meta analysis of the top features added
## one at a time
library(ggplot2)
## let's load the results of the meta-analysis
load('WO_meta_26C2.rda')
col=c('firebrick2','darkcyan','forestgreen','darkorange1')
## Now let's plot this

xx2 <- 1:20
nSubj   <- 26                  # how many subjects in your dataset?
nFeats  <- 20                  # max number of features on x-axis of final plot
x       <- 2                   # how many do you want to leave out?
Folds   <- combn(nSubj,x)      # all the possible combination of Leave x out
nReps   <- dim(Folds)[2]       # make sure to cover all possibilities

# size of matrices is nReps x nFeats

xx <- t(kronecker(matrix(1,1,nReps),1:nFeats))
plot(xx[,1],WO_LS_L[,1], ylim = c(0, 1), xlim = c(0.5, nFeats+0.5), col='firebrick2',
     xlab=paste('Features'),ylab=paste('R Squared'),cex.lab=0.8,pch=21, cex=0.8,lwd=2)
points(xx[,1],WO_LS_Q[,1], col = 'darkcyan',pch=22, cex=0.8)

for (i in 2:nFeats)
{
  points(xx[,i],WO_LS_L[,i], col='firebrick2',pch=21, cex=0.8)
  points(xx[,i],WO_LS_Q[,i], col = 'darkcyan',pch=22, cex=0.8)
  
}
yy <- colMeans(WO_LS_L)
yy2 <- colMeans(WO_LS_Q)

## now let's plot smooth curves through the mean R^2 values:
lo <- loess(yy~xx2)
lines(predict(lo), col='firebrick2', lwd=2)
lo <- loess(yy2~xx2)
lines(predict(lo), col='darkcyan', lwd=2)
