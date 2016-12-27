% This function compares the RMSE and R2 values for LASSO analyses
% performend using different methods (see version control file for details
% (opens in both matlab and R) on the methods used. This will address 
% v8-v13 in LASSO and v8-v10 in random forests
% This is only for RMSE, for the other metrics, they have to be summarized
% and are thus in the compute_performance_metrics function

clear all
clc
%% Load all the files
% the models are in order: Lin, Quad with no extrapolation
% load all model results (all predictions are in the form nSubj x nReps
% nSubj: Number of Subjects
% nReps: Number of cross-validation repeats
% change in this, dec 14, 2016:
% edited to reflect new csv names, removed references to weak-performing
% lasso (externally validated lasso, v8-v10). Only use best lasso can do
% (v11-v13)
% v8, v11: nReps 100, nFolds 4
% v9, v12: nReps 1, nFolds 26 (LOO)
% v10, v13: nReps 100, nFolds 13

FM_RF_L_v8 = loadmodres('2016dec12_FM_RF_L_nReps100_nFolds4.csv'); 
FM_RF_Q_v8 = loadmodres('2016dec12_FM_RF_Q_nReps100_nFolds4.csv');

WO_RF_L_v8 = loadmodres('2016dec12_WO_RF_L_nReps100_nFolds4.csv'); 
WO_RF_Q_v8 = loadmodres('2016dec12_WO_RF_Q_nReps100_nFolds4.csv');

FM_RF_L_v9 = loadmodres('2016dec12_FM_RF_L_nReps1_nFolds26.csv'); 
FM_RF_Q_v9 = loadmodres('2016dec12_FM_RF_Q_nReps1_nFolds26.csv');

WO_RF_L_v9 = loadmodres('2016dec12_WO_RF_L_nReps1_nFolds26.csv'); 
WO_RF_Q_v9 = loadmodres('2016dec12_WO_RF_Q_nReps1_nFolds26.csv');

FM_RF_L_v10 = loadmodres('2016dec12_FM_RF_L_nReps100_nFolds13.csv'); 
FM_RF_Q_v10 = loadmodres('2016dec12_FM_RF_Q_nReps100_nFolds13.csv');

WO_RF_L_v10 = loadmodres('2016dec12_WO_RF_L_nReps100_nFolds13.csv'); 
WO_RF_Q_v10 = loadmodres('2016dec12_WO_RF_Q_nReps100_nFolds13.csv');

FM_LS_L_v11 = loadmodres('2016dec12_FM_LS_L_nReps100_nFolds4.csv'); 
FM_LS_Q_v11 = loadmodres('2016dec12_FM_LS_Q_nReps100_nFolds4.csv'); 

WO_LS_L_v11 = loadmodres('2016dec12_WO_LS_L_nReps100_nFolds4.csv'); 
WO_LS_Q_v11 = loadmodres('2016dec12_WO_LS_Q_nReps100_nFolds4.csv');

FM_LS_L_v12 = loadmodres('2016dec12_FM_LS_L_nReps1_nFolds26.csv'); 
FM_LS_Q_v12 = loadmodres('2016dec12_FM_LS_Q_nReps1_nFolds26.csv'); 

WO_LS_L_v12 = loadmodres('2016dec12_WO_LS_L_nReps1_nFolds26.csv'); 
WO_LS_Q_v12 = loadmodres('2016dec12_WO_LS_Q_nReps1_nFolds26.csv'); 

FM_LS_L_v13 = loadmodres('2016dec12_FM_LS_L_nReps100_nFolds13.csv'); 
FM_LS_Q_v13 = loadmodres('2016dec12_FM_LS_Q_nReps100_nFolds13.csv'); 

WO_LS_L_v13 = loadmodres('2016dec12_WO_LS_L_nReps100_nFolds13.csv'); 
WO_LS_Q_v13 = loadmodres('2016dec12_WO_LS_Q_nReps100_nFolds13.csv'); 

fid = fopen('yFM.csv'); pit = textscan(fid,'%s%s','delimiter',','); fclose(fid); yFM = str2double(pit{2}(2:end));
fid = fopen('yWO.csv'); pit = textscan(fid,'%s%s','delimiter',','); fclose(fid); yWO = str2double(pit{2}(2:end));

clearvars pit fid ans

% Remember that v9 and v12 are LOO-CV, so they have only 1 repeat. This
% means the matrices need to be reduced to one column

FM_RF_L_v9(:,2:end) = []; FM_RF_Q_v9(:,2:end) = [];
WO_RF_L_v9(:,2:end) = []; WO_RF_Q_v9(:,2:end) = [];
WO_LS_L_v12(:,2:end) = []; WO_LS_Q_v12(:,2:end) = [];

%% Intialize variables to hold RMSE and R2
nReps = size(FM_RF_L_v8,2); % number of repeats of CV for 4 and 13-fold models, LOO models have a nReps of 1

% Will now construct RMSE, R^2, adjusted R^2 and Slope as cell arrays, with
% the following order, one column for each FM and WO, first column is FM
% 01- LS_L_v8 
% 02- LS_Q_v8
% 03- LS_L_v9
% 04- LS_Q_v9
% 05- LS_L_v10
% 06- LS_Q_v10
% 07- LS_L_v11
% 08- LS_Q_v11
% 09- LS_L_v12
% 10- LS_Q_v12
% 11- LS_L_v13
% 12- LS_Q_v13
% 13- RF_L_v8 
% 14- RF_Q_v8
% 15- RF_L_v9
% 16- RF_Q_v9
% 17- RF_L_v10
% 18- RF_Q_v10
RMSE = cell(18,2);

%% find RMSE values for each subject for each cv repeat

for i = 1:100 % loop over columns
    RMSE{7,1}(:,i) = abs(FM_LS_L_v11(:,i) - yFM);       RMSE{7,2}(:,i) = abs(WO_LS_L_v11(:,i) - yWO);
    RMSE{8,1}(:,i) = abs(FM_LS_Q_v11(:,i) - yFM);       RMSE{8,2}(:,i) = abs(WO_LS_Q_v11(:,i) - yWO);
    RMSE{11,1}(:,i) = abs(FM_LS_L_v13(:,i) - yFM);      RMSE{11,2}(:,i) = abs(WO_LS_L_v13(:,i) - yWO);
    RMSE{12,1}(:,i) = abs(FM_LS_Q_v13(:,i) - yFM);      RMSE{12,2}(:,i) = abs(WO_LS_Q_v13(:,i) - yWO);
    RMSE{13,1}(:,i) = abs(FM_RF_L_v8(:,i) - yFM);       RMSE{13,2}(:,i) = abs(WO_RF_L_v8(:,i) - yWO);
    RMSE{14,1}(:,i) = abs(FM_RF_Q_v8(:,i) - yFM);       RMSE{14,2}(:,i) = abs(WO_RF_Q_v8(:,i) - yWO);
    RMSE{17,1}(:,i) = abs(FM_RF_L_v10(:,i) - yFM);      RMSE{17,2}(:,i) = abs(WO_RF_L_v10(:,i) - yWO);
    RMSE{18,1}(:,i) = abs(FM_RF_Q_v10(:,i) - yFM);      RMSE{18,2}(:,i) = abs(WO_RF_Q_v10(:,i) - yWO);
end
RMSE{9,1}(:,1) = abs(FM_LS_L_v12(:,1) - yFM);   RMSE{9,2}(:,1) = abs(WO_LS_L_v12(:,1) - yWO);
RMSE{10,1}(:,1) = abs(FM_LS_Q_v12(:,1) - yFM);  RMSE{10,2}(:,1) = abs(WO_LS_Q_v12(:,1) - yWO);
RMSE{15,1}(:,1) = abs(FM_RF_L_v9(:,1) - yFM);   RMSE{15,2}(:,1) = abs(WO_RF_L_v9(:,1) - yWO);
RMSE{16,1}(:,1) = abs(FM_RF_Q_v9(:,1) - yFM);   RMSE{16,2}(:,1) = abs(WO_RF_Q_v9(:,1) - yWO);

% clearvars -except RMSE R2 AR2 SLP yFM yWO

%% Plot v8 (RF) and v11(LASSO)
col= zeros(4,3);
col(1,:) = [1 0 0]; col(2,:) = [0 0 1]; col(3,:) = [1 0.7 0]; col(4,:) = [0 1 0];

figure
clf
axes('Position',[0.05 0.02 0.95 0.95])
hold on

scf = 1;
plot1DDistributionV2(RMSE{7,1}(:),'w',[0.8,1.2],col(1,:),scf)
plot1DDistributionV2(RMSE{8,1}(:),'w',[1.8,2.2],col(2,:),scf)
plot1DDistributionV2(RMSE{13,1}(:),'w',[2.8,3.2],col(3,:),scf)
plot1DDistributionV2(RMSE{14,1}(:),'w',[3.8,4.2],col(4,:),scf)

plot1DDistributionV2(RMSE{7,2}(:),'w',[5.8,6.5],col(1,:),scf)
plot1DDistributionV2(RMSE{8,2}(:),'w',[6.8,7.2],col(2,:),scf)
plot1DDistributionV2(RMSE{13,2}(:),'w',[7.8,8.2],col(3,:),scf)
plot1DDistributionV2(RMSE{14,2}(:),'w',[8.8,9.2],col(4,:),scf)

ax = gca;
% ax.XTick = [2.5 7.5];
ax.XTick = [];
% ax.XTickLabel = {'Fugl-Meyer','Wolf Motor Function'};
set(ax,'FontSize',18)
ylabel('RMSE')
ylim([0 28])
xlim([0 10])
