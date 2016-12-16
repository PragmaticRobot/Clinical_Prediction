
# This function extracts the list of features that were given a coefficient by the lasso
# algorithm, these form the shortlist of features needed to make a prediction


#### V8 ####

# load('LS_v8.rda')
# 
# nReps = length(FM_LS_L_beta)
# nFolds = length(FM_LS_L_beta[[1]])
# nFeat = length(FM_LS_L_beta[[1]][[1]])
# 
# FM_LS_L_coef_v8 <- as.data.frame(matrix(0,nFeat,nReps*nFolds))
# rownames(FM_LS_L_coef_v8) <- names(FM_LS_L_beta[[1]][[1]])
# 
# for (i in 1:nReps)
# {
#   for (kk in 1:nFolds)
#   {
#     coln <- (nFolds*(i-1))+kk # column number
#     FM_LS_L_coef_v8[,coln] <- FM_LS_L_beta[[i]][[kk]]
#   }
# }
# 
# counts_v8 <- rowSums(FM_LS_L_coef_v8 != 0)
# counts_v8 <- as.data.frame(counts_v8)

#### V9 ####
# 
# load('LS_v9.rda')
# 
# nReps = length(FM_LS_L_beta)
# nFolds = length(FM_LS_L_beta[[1]])
# nFeat = length(FM_LS_L_beta[[1]][[1]])
# 
# FM_LS_L_coef_v9 <- as.data.frame(matrix(0,nFeat,nReps*nFolds))
# rownames(FM_LS_L_coef_v9) <- names(FM_LS_L_beta[[1]][[1]])
# 
# for (i in 1:nReps)
# {
#   for (kk in 1:nFolds)
#   {
#     coln <- (nFolds*(i-1))+kk # column number
#     FM_LS_L_coef_v9[,coln] <- FM_LS_L_beta[[i]][[kk]]
#   }
# }
# 
# counts_v9 <- rowSums(FM_LS_L_coef_v9 != 0)
# counts_v9 <- as.data.frame(counts_v9)

#### V10 ####
# 
# load('LS_v10.rda')
# 
# nReps = length(FM_LS_L_beta)
# nFolds = length(FM_LS_L_beta[[1]])
# nFeat = length(FM_LS_L_beta[[1]][[1]])
# 
# FM_LS_L_coef_v10 <- as.data.frame(matrix(0,nFeat,nReps*nFolds))
# rownames(FM_LS_L_coef_v10) <- names(FM_LS_L_beta[[1]][[1]])
# 
# for (i in 1:nReps)
# {
#   for (kk in 1:nFolds)
#   {
#     coln <- (nFolds*(i-1))+kk # column number
#     FM_LS_L_coef_v10[,coln] <- FM_LS_L_beta[[i]][[kk]]
#   }
# }
# 
# counts_v10 <- rowSums(FM_LS_L_coef_v10 != 0)
# counts_v10 <- as.data.frame(counts_v10)
#### v11 ####
source('create_df1_df2_no_unknown.R')
load('2016dec12_LS_nReps100_nFolds4.rda')

nReps = length(FM_LS_L_beta)
nFolds = 1 # even though this is 4 fold cv, for each repeat only 1 model is used
nFeat = 52 # 52 features + intercept

FM_LS_L_coef_v11 <- as.data.frame(matrix(0,nFeat,nReps))
rownames(FM_LS_L_coef_v11) <- c(gverb2)

for (jj in 1:nReps)
{
  FM_LS_L_coef_v11[FM_LS_L_beta[[jj]]$i,jj] <- FM_LS_L_beta[[jj]]$x
}

counts_v11 <- rowSums(FM_LS_L_coef_v11 != 0)
counts_v11 <- as.data.frame(counts_v11)

#### v12 ####

load('2016dec12_LS_nReps1_nFolds26.rda')

nReps = length(FM_LS_L_beta)

FM_LS_L_coef_v12 <- as.data.frame(matrix(0,nFeat,nReps))
rownames(FM_LS_L_coef_v12) <- c(gverb2)

for (jj in 1:nReps)
{
  FM_LS_L_coef_v12[FM_LS_L_beta[[jj]]$i,jj] <- FM_LS_L_beta[[jj]]$x
}

counts_v12 <- rowSums(FM_LS_L_coef_v12 != 0)
counts_v12 <- as.data.frame(counts_v12)
#### v13 ####

load('2016dec12_LS_nReps100_nFolds13.rda')

nReps = length(FM_LS_L_beta)

FM_LS_L_coef_v13 <- as.data.frame(matrix(0,nFeat,nReps))
rownames(FM_LS_L_coef_v13) <- c(gverb2)

for (jj in 1:nReps)
{
  FM_LS_L_coef_v13[FM_LS_L_beta[[jj]]$i,jj] <- FM_LS_L_beta[[jj]]$x
}

counts_v13 <- rowSums(FM_LS_L_coef_v13 != 0)
counts_v13 <- as.data.frame(counts_v13)

# next remove the first row (intercept) from all counts
counts_v11 <- counts_v11[-c(1), , drop=FALSE]
counts_v12 <- counts_v12[-c(1), , drop=FALSE]
counts_v13 <- counts_v13[-c(1), , drop=FALSE]

# counts_total <- cbind(counts_v8,counts_v9,counts_v10, counts_v11, counts_v12, counts_v13)
counts_total <- cbind(counts_v11, counts_v12, counts_v13)
# now we need to normalize counts so we can compare

# counts_total$counts_v8 <- counts_total$counts_v8/400
# counts_total$counts_v9 <- counts_total$counts_v9/26
# counts_total$counts_v10 <- counts_total$counts_v10/1300
counts_total$counts_v11 <- counts_total$counts_v11/100
counts_total$counts_v12 <- counts_total$counts_v12/1
counts_total$counts_v13 <- counts_total$counts_v13/100

write.csv(counts_total,"2016dec12_FM_lasso_counts.csv")

# save(counts_v8,counts_v9,counts_v10,counts_v11,counts_v12,counts_v13,counts_total,
     # file = 'lasso_FM_counts.rda')
save(counts_v11,counts_v12,counts_v13,counts_total,
     file = '2016dec12_lasso_FM_counts.rda')