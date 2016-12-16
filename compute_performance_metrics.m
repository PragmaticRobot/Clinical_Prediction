% This function compares the RMSE and R2 values for LASSO analyses
% performend using different methods (see version control file for details
% (opens in both matlab and R) on the methods used. This will address 
% v8-v13 in LASSO and v8-v10 in random forests

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
R2 = cell(18,2);
AR2 = cell(18,2);
SLP = cell(18,2);

%% Get all RMSE and Rsquared and slopes
for i = 1:nReps
    trash = fitlm(FM_LS_L_v11(:,i),yFM); RMSE{7,1}(i) = trash.RMSE; R2{7,1}(i) = trash.Rsquared.Ordinary; AR2{7,1}(i) = trash.Rsquared.Adjusted; SLP{7,1}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_LS_Q_v11(:,i),yFM); RMSE{8,1}(i) = trash.RMSE; R2{8,1}(i) = trash.Rsquared.Ordinary; AR2{8,1}(i) = trash.Rsquared.Adjusted; SLP{8,1}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_LS_L_v13(:,i),yFM); RMSE{11,1}(i) = trash.RMSE; R2{11,1}(i) = trash.Rsquared.Ordinary; AR2{11,1}(i) = trash.Rsquared.Adjusted; SLP{11,1}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_LS_Q_v13(:,i),yFM); RMSE{12,1}(i) = trash.RMSE; R2{12,1}(i) = trash.Rsquared.Ordinary; AR2{12,1}(i) = trash.Rsquared.Adjusted; SLP{12,1}(i) = trash.Coefficients.Estimate(2);
    
    trash = fitlm(FM_RF_L_v8(:,i),yFM); RMSE{13,1}(i) = trash.RMSE; R2{13,1}(i) = trash.Rsquared.Ordinary; AR2{13,1}(i) = trash.Rsquared.Adjusted; SLP{13,1}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_RF_Q_v8(:,i),yFM); RMSE{14,1}(i) = trash.RMSE; R2{14,1}(i) = trash.Rsquared.Ordinary; AR2{14,1}(i) = trash.Rsquared.Adjusted; SLP{14,1}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_RF_L_v10(:,i),yFM); RMSE{17,1}(i) = trash.RMSE; R2{17,1}(i) = trash.Rsquared.Ordinary; AR2{17,1}(i) = trash.Rsquared.Adjusted; SLP{17,1}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_RF_Q_v10(:,i),yFM); RMSE{18,1}(i) = trash.RMSE; R2{18,1}(i) = trash.Rsquared.Ordinary; AR2{18,1}(i) = trash.Rsquared.Adjusted; SLP{18,1}(i) = trash.Coefficients.Estimate(2);
    
    trash = fitlm(WO_LS_L_v11(:,i),yWO); RMSE{7,2}(i) = trash.RMSE; R2{7,2}(i) = trash.Rsquared.Ordinary; AR2{7,2}(i) = trash.Rsquared.Adjusted; SLP{7,2}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_LS_Q_v11(:,i),yWO); RMSE{8,2}(i) = trash.RMSE; R2{8,2}(i) = trash.Rsquared.Ordinary; AR2{8,2}(i) = trash.Rsquared.Adjusted; SLP{8,2}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_LS_L_v13(:,i),yWO); RMSE{11,2}(i) = trash.RMSE; R2{11,2}(i) = trash.Rsquared.Ordinary; AR2{11,2}(i) = trash.Rsquared.Adjusted; SLP{11,2}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_LS_Q_v13(:,i),yWO); RMSE{12,2}(i) = trash.RMSE; R2{12,2}(i) = trash.Rsquared.Ordinary; AR2{12,2}(i) = trash.Rsquared.Adjusted; SLP{12,2}(i) = trash.Coefficients.Estimate(2);
    
    trash = fitlm(WO_RF_L_v8(:,i),yWO); RMSE{13,2}(i) = trash.RMSE; R2{13,2}(i) = trash.Rsquared.Ordinary; AR2{13,2}(i) = trash.Rsquared.Adjusted; SLP{13,2}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_RF_Q_v8(:,i),yWO); RMSE{14,2}(i) = trash.RMSE; R2{14,2}(i) = trash.Rsquared.Ordinary; AR2{14,2}(i) = trash.Rsquared.Adjusted; SLP{14,2}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_RF_L_v10(:,i),yWO); RMSE{17,2}(i) = trash.RMSE; R2{17,2}(i) = trash.Rsquared.Ordinary; AR2{17,2}(i) = trash.Rsquared.Adjusted; SLP{17,2}(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_RF_Q_v10(:,i),yWO); RMSE{18,2}(i) = trash.RMSE; R2{18,2}(i) = trash.Rsquared.Ordinary; AR2{18,2}(i) = trash.Rsquared.Adjusted; SLP{18,2}(i) = trash.Coefficients.Estimate(2);
end

trash = fitlm(FM_LS_L_v9(:,1),yFM); RMSE{3,1}(1) = trash.RMSE; R2{3,1}(1) = trash.Rsquared.Ordinary; AR2{3,1}(1) = trash.Rsquared.Adjusted; SLP{3,1}(1) = trash.Coefficients.Estimate(2);
trash = fitlm(FM_LS_Q_v9(:,1),yFM); RMSE{4,1}(1) = trash.RMSE; R2{4,1}(1) = trash.Rsquared.Ordinary; AR2{4,1}(1) = trash.Rsquared.Adjusted; SLP{4,1}(1) = trash.Coefficients.Estimate(2);

trash = fitlm(FM_LS_L_v12(:,1),yFM); RMSE{9,1}(1) = trash.RMSE; R2{9,1}(1) = trash.Rsquared.Ordinary; AR2{9,1}(1) = trash.Rsquared.Adjusted; SLP{9,1}(1) = trash.Coefficients.Estimate(2);
trash = fitlm(FM_LS_Q_v12(:,1),yFM); RMSE{10,1}(1) = trash.RMSE; R2{10,1}(1) = trash.Rsquared.Ordinary; AR2{10,1}(1) = trash.Rsquared.Adjusted; SLP{10,1}(1) = trash.Coefficients.Estimate(2);

trash = fitlm(FM_RF_L_v9(:,1),yFM); RMSE{15,1}(1) = trash.RMSE; R2{15,1}(1) = trash.Rsquared.Ordinary; AR2{15,1}(1) = trash.Rsquared.Adjusted; SLP{15,1}(1) = trash.Coefficients.Estimate(2);
trash = fitlm(FM_RF_Q_v9(:,1),yFM); RMSE{16,1}(1) = trash.RMSE; R2{16,1}(1) = trash.Rsquared.Ordinary; AR2{16,1}(1) = trash.Rsquared.Adjusted; SLP{16,1}(1) = trash.Coefficients.Estimate(2);

trash = fitlm(WO_LS_L_v9(:,1),yWO); RMSE{3,2}(1) = trash.RMSE; R2{3,2}(1) = trash.Rsquared.Ordinary; AR2{3,2}(1) = trash.Rsquared.Adjusted; SLP{3,2}(1) = trash.Coefficients.Estimate(2);
trash = fitlm(WO_LS_Q_v9(:,1),yWO); RMSE{4,2}(1) = trash.RMSE; R2{4,2}(1) = trash.Rsquared.Ordinary; AR2{4,2}(1) = trash.Rsquared.Adjusted; SLP{4,2}(1) = trash.Coefficients.Estimate(2);
    
trash = fitlm(WO_LS_L_v12(:,1),yWO); RMSE{9,2}(1) = trash.RMSE; R2{9,2}(1) = trash.Rsquared.Ordinary; AR2{9,2}(1) = trash.Rsquared.Adjusted; SLP{9,2}(1) = trash.Coefficients.Estimate(2);
trash = fitlm(WO_LS_Q_v12(:,1),yWO); RMSE{10,2}(1) = trash.RMSE; R2{10,2}(1) = trash.Rsquared.Ordinary; AR2{10,2}(1) = trash.Rsquared.Adjusted; SLP{10,2}(1) = trash.Coefficients.Estimate(2);

trash = fitlm(WO_RF_L_v9(:,1),yWO); RMSE{15,2}(1) = trash.RMSE; R2{15,2}(1) = trash.Rsquared.Ordinary; AR2{15,2}(1) = trash.Rsquared.Adjusted; SLP{15,2}(1) = trash.Coefficients.Estimate(2);
trash = fitlm(WO_RF_Q_v9(:,1),yWO); RMSE{16,2}(1) = trash.RMSE; R2{16,2}(1) = trash.Rsquared.Ordinary; AR2{16,2}(1) = trash.Rsquared.Adjusted; SLP{16,2}(1) = trash.Coefficients.Estimate(2);

clearvars -except RMSE R2 AR2 SLP yFM yWO

%% calculate mean and std for all model types

RMSE_M = zeros(18,2); RMSE_STD = zeros(18,2); R2_M = zeros(18,2); R2_STD = zeros(18,2);
AR2_M = zeros(18,2); AR2_STD = zeros(18,2); SLP_M = zeros(18,2); SLP_STD = zeros(18,2);

for i = 1:2        % loop over columns
    for j = 1:18        % loop over rows
        RMSE_M(j,i) = mean(RMSE{j,i});
        R2_M(j,i) = mean(R2{j,i});
        AR2_M(j,i) = mean(AR2{j,i});
        SLP_M(j,i) = mean(SLP{j,i});
        
        RMSE_STD(j,i) = std(RMSE{j,i});
        R2_STD(j,i) = std(R2{j,i});
        AR2_STD(j,i) = std(AR2{j,i});
        SLP_STD(j,i) = std(SLP{j,i});
    end
end

save('2016dec12_performance_metrics.mat')

%% Let's get the bar graphs
% plot_comparison_bars(RMSE_M,RMSE_STD,'RMSE')
% plot_comparison_bars(R2_M,R2_STD,'R2')
% plot_comparison_bars(AR2_M,AR2_STD,'Adjusted R2')
% plot_comparison_bars(SLP_M,SLP_STD,'Slope')
