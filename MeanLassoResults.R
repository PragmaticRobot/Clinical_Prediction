load('lasso_FM_counts.rda')

FM_counts <- counts_total

load('lasso_WO_counts.rda')

WO_counts <- counts_total

FM_tots <- data.frame(Means = rowMeans(FM_counts))
WO_tots <- data.frame(Means = rowMeans(WO_counts))

write.csv(df1, file = 'df1.csv')
write.csv(df3, file = 'df3.csv')
write.csv(yFM, file = 'yFM.csv')
write.csv(yWO, file = 'yWO.csv')

