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

FM_LS_L_v8 = loadmodres('mod1_v8.csv'); FM_LS_Q_v8 = loadmodres('mod2_v8.csv'); FM_RF_L_v8 = loadmodres('mod3_v8.csv'); FM_RF_Q_v8 = loadmodres('mod4_v8.csv');
WO_LS_L_v8 = loadmodres('mod5_v8.csv'); WO_LS_Q_v8 = loadmodres('mod6_v8.csv'); WO_RF_L_v8 = loadmodres('mod7_v8.csv'); WO_RF_Q_v8 = loadmodres('mod8_v8.csv');

FM_LS_L_v9 = loadmodres('mod1_v9.csv'); FM_LS_Q_v9 = loadmodres('mod2_v9.csv'); FM_RF_L_v9 = loadmodres('mod3_v9.csv'); FM_RF_Q_v9 = loadmodres('mod4_v9.csv');
WO_LS_L_v9 = loadmodres('mod5_v9.csv'); WO_LS_Q_v9 = loadmodres('mod6_v9.csv'); WO_RF_L_v9 = loadmodres('mod7_v9.csv'); WO_RF_Q_v9 = loadmodres('mod8_v9.csv');

FM_LS_L_v10 = loadmodres('mod1_v10.csv'); FM_LS_Q_v10 = loadmodres('mod2_v10.csv'); FM_RF_L_v10 = loadmodres('mod3_v10.csv'); FM_RF_Q_v10 = loadmodres('mod4_v10.csv');
WO_LS_L_v10 = loadmodres('mod5_v10.csv'); WO_LS_Q_v10 = loadmodres('mod6_v10.csv'); WO_RF_L_v10 = loadmodres('mod7_v10.csv'); WO_RF_Q_v10 = loadmodres('mod8_v10.csv');

FM_LS_L_v11 = loadmodres('mod1_v11.csv'); FM_LS_Q_v11 = loadmodres('mod2_v11.csv'); WO_LS_L_v11 = loadmodres('mod5_v11.csv'); WO_LS_Q_v11 = loadmodres('mod6_v11.csv');

FM_LS_L_v12 = loadmodres('mod1_v12.csv'); FM_LS_Q_v12 = loadmodres('mod2_v12.csv'); WO_LS_L_v12 = loadmodres('mod5_v12.csv'); WO_LS_Q_v12 = loadmodres('mod6_v12.csv'); 

FM_LS_L_v13 = loadmodres('mod1_v13.csv'); FM_LS_Q_v13 = loadmodres('mod2_v13.csv'); WO_LS_L_v13 = loadmodres('mod5_v13.csv'); WO_LS_Q_v13 = loadmodres('mod6_v13.csv'); 

fid = fopen('yFM.csv'); pit = textscan(fid,'%s%s','delimiter',','); fclose(fid); yFM = str2double(pit{2}(2:end));
fid = fopen('yWO.csv'); pit = textscan(fid,'%s%s','delimiter',','); fclose(fid); yWO = str2double(pit{2}(2:end));

clearvars pit fid ans

% Remember that v9 and v12 are LOO-CV, so they have only 1 repeat. This
% means the matrices need to be reduced to one column

FM_LS_L_v9(:,2:end) = []; FM_LS_Q_v9(:,2:end) = []; FM_RF_L_v9(:,2:end) = []; FM_RF_Q_v9(:,2:end) = [];
WO_LS_L_v9(:,2:end) = []; WO_LS_Q_v9(:,2:end) = []; WO_RF_L_v9(:,2:end) = []; WO_RF_Q_v9(:,2:end) = [];
FM_LS_L_v12(:,2:end) = []; FM_LS_Q_v12(:,2:end) = []; WO_LS_L_v12(:,2:end) = []; WO_LS_Q_v12(:,2:end) = [];

%% Intialize variables to hold RMSE and R2
nReps = size(FM_LS_L_v8,2); % number of repeats of CV for 4 and 13-fold models, LOO models have a nReps of 1

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
    RMSE{1,1}(:,i) = abs(FM_LS_L_v8(:,i) - yFM);        RMSE{1,2}(:,i) = abs(WO_LS_L_v8(:,i) - yFM);
    RMSE{2,1}(:,i) = abs(FM_LS_Q_v8(:,i) - yFM);        RMSE{2,2}(:,i) = abs(WO_LS_Q_v8(:,i) - yFM);
    RMSE{5,1}(:,i) = abs(FM_LS_L_v10(:,i) - yFM);       RMSE{5,2}(:,i) = abs(WO_LS_L_v10(:,i) - yFM);
    RMSE{6,1}(:,i) = abs(FM_LS_Q_v10(:,i) - yFM);       RMSE{6,2}(:,i) = abs(WO_LS_Q_v10(:,i) - yFM);
    RMSE{7,1}(:,i) = abs(FM_LS_L_v11(:,i) - yFM);       RMSE{7,2}(:,i) = abs(WO_LS_L_v11(:,i) - yFM);
    RMSE{8,1}(:,i) = abs(FM_LS_Q_v11(:,i) - yFM);       RMSE{8,2}(:,i) = abs(WO_LS_Q_v11(:,i) - yFM);
    RMSE{11,1}(:,i) = abs(FM_LS_L_v13(:,i) - yFM);      RMSE{11,2}(:,i) = abs(WO_LS_L_v13(:,i) - yFM);
    RMSE{12,1}(:,i) = abs(FM_LS_Q_v13(:,i) - yFM);      RMSE{12,2}(:,i) = abs(WO_LS_Q_v13(:,i) - yFM);
    RMSE{13,1}(:,i) = abs(FM_RF_L_v8(:,i) - yFM);       RMSE{13,2}(:,i) = abs(WO_RF_L_v8(:,i) - yFM);
    RMSE{14,1}(:,i) = abs(FM_RF_Q_v8(:,i) - yFM);       RMSE{14,2}(:,i) = abs(WO_RF_Q_v8(:,i) - yFM);
    RMSE{17,1}(:,i) = abs(FM_RF_L_v10(:,i) - yFM);      RMSE{17,2}(:,i) = abs(WO_RF_L_v10(:,i) - yFM);
    RMSE{18,1}(:,i) = abs(FM_RF_Q_v10(:,i) - yFM);      RMSE{18,2}(:,i) = abs(WO_RF_Q_v10(:,i) - yFM);
end
RMSE{3,1}(:,1) = abs(FM_LS_L_v9(:,1) - yFM);    RMSE{3,2}(:,1) = abs(WO_LS_L_v9(:,1) - yFM);
RMSE{4,1}(:,1) = abs(FM_LS_Q_v9(:,1) - yFM);    RMSE{4,2}(:,1) = abs(WO_LS_Q_v9(:,1) - yFM);
RMSE{9,1}(:,1) = abs(FM_LS_L_v12(:,1) - yFM);   RMSE{9,2}(:,1) = abs(WO_LS_L_v12(:,1) - yFM);
RMSE{10,1}(:,1) = abs(FM_LS_Q_v12(:,1) - yFM);  RMSE{10,2}(:,1) = abs(WO_LS_Q_v12(:,1) - yFM);
RMSE{15,1}(:,1) = abs(FM_RF_L_v9(:,1) - yFM);   RMSE{15,2}(:,1) = abs(WO_RF_L_v9(:,1) - yFM);
RMSE{16,1}(:,1) = abs(FM_RF_Q_v9(:,1) - yFM);   RMSE{16,2}(:,1) = abs(WO_RF_Q_v9(:,1) - yFM);

% clearvars -except RMSE R2 AR2 SLP yFM yWO

%% calculate mean and std for all model types

RMSE_M = zeros(18,2); RMSE_STD = zeros(18,2); R2_M = zeros(18,2); R2_STD = zeros(18,2);
AR2_M = zeros(18,2); AR2_STD = zeros(18,2); SLP_M = zeros(18,2); SLP_STD = zeros(18,2);

for i = 1:2        % loop over columns
    for j = 1:18        % loop over rows
        This_RMSE = RMSE{j,i}(:);
        RMSE_M(j,i) = mean(This_RMSE);
        RMSE_STD(j,i) = std(This_RMSE);
    end
end

save('UnsummarizedRMSE.mat')
% 
% locs = [1,2,4,5,7,8,12,13,15,16,18,19];
% means = RMSE_M(1:12,1);
% stds = RMSE_STD(1:12,1);
% plot_performance_bars(means, stds,locs,1,1,1)
% 
% % Now let's compare best lasso and RF performance in both RMSE and R2
% locs = [1,2,3,4,6,7,8,9,11,12,13,14];
% means = RMSE_M([7,8,13,14,9,10,15,16,11,12,17,18],1);
% stds = RMSE_STD([7,8,13,14,9,10,15,16,11,12,17,18],1);
% plot_performance_bars(means, stds,locs,1,2,1)
% 
% % Do the same for Wolf
% locs = [1,2,4,5,7,8,12,13,15,16,18,19];
% means = RMSE_M(1:12,2);
% stds = RMSE_STD(1:12,2);
% plot_performance_bars(means, stds,locs,1,1,2)
% 
% % Now let's compare best lasso and RF performance in both RMSE and R2
% locs = [1,2,3,4,6,7,8,9,11,12,13,14];
% means = RMSE_M([7,8,13,14,9,10,15,16,11,12,17,18],2);
% stds = RMSE_STD([7,8,13,14,9,10,15,16,11,12,17,18],2);
% plot_performance_bars(means, stds,locs,1,2,2)
figure
% hold on
% for i = 1:18
%     for j = 1:2
%         subplot(9,4,i); hold on
%         plot1DDistribution(RMSE{i,j}(:),'b')
%     end
% end
subplot(1,2,1); plot1DDistribution(RMSE{9,1}(:))
subplot(1,2,2); plot1DDistribution(RMSE{7,1}(:))
