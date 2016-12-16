
## This script calcuates the mean lasso results, and loads the RF VI results.
## Then 

rm(list = ls()) #clear the workspace
load('2016dec12_Basics.rda')

## load and clean up lasso
load('2016dec12_lasso_FM_counts.rda')
FM_counts <- counts_total

load('2016dec12_lasso_WO_counts.rda')
WO_counts <- counts_total

LS_FM_tots <- data.frame(Means = rowMeans(FM_counts))
LS_WO_tots <- data.frame(Means = rowMeans(WO_counts))

## load and clean up RF

load('RF_Ranks_nReps100_nFolds4.rda')
RF_FM_v10 <- ranks_FM
RF_WO_v10 <- ranks_WO

load('RF_Ranks_nReps1_nFolds26.rda')
RF_FM_v9 <- ranks_FM
RF_WO_v9 <- ranks_WO

load('RF_Ranks_nReps100_nFolds13.rda')
RF_FM_v8 <- ranks_FM
RF_WO_v8 <- ranks_WO

rm(counts_v8, counts_v9, counts_v10, counts_v11, counts_v12, counts_v13, counts_total,
   VIs_FM,VIs_WO, ranks_FM,ranks_WO, FM_counts, WO_counts)

write.csv(LS_FM_tots,file = '2016dec12_LS_FM_order.csv')
write.csv(LS_WO_tots,file = '2016dec12_LS_WO_order.csv')
write.csv(RF_FM_v9,file = 'RF_FM_order.csv')
write.csv(RF_WO_v9,file = 'RF_WO_order.csv')
