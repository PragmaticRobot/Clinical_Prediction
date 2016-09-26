clear all
clc
%% Load all the files
% the models are in order: Lin, NoExtrap, Quad
% 1-3: FM LS
% 4-6: FM RF
% 7-9: WO LS
% 10-12: WO RF
FM_Lin_LS = loadmodres('FM_lin_LASSO_4FoldNoExtrap.csv');
FM_NoExtrap_LS = loadmodres('FM_NoExtrap_LASSO_4FoldNoExtrap.csv');
FM_quad_LS = loadmodres('FM_quad_LASSO_4FoldNoExtrap.csv');
FM_Lin_RF = loadmodres('FM_Lin_RF_4FoldNoExtrap.csv');
FM_NoExtrap_RF = loadmodres('FM_NoExtrap_RF_4FoldNoExtrap.csv');
FM_quad_RF = loadmodres('FM_quad_RF_4FoldNoExtrap.csv');
