clear all
clc
%% Load all the files
% the models are in order: Lin, NoExtrap, Quad
% 1-2: FM LS
% 4-5: FM RF
% 7-8: WO LS
% 10-11: WO RF
FM_RF_L = loadmodres('2016dec12_FM_RF_L_nReps100_nFolds4.csv'); 
FM_RF_Q = loadmodres('2016dec12_FM_RF_Q_nReps100_nFolds4.csv');

WO_RF_L = loadmodres('2016dec12_WO_RF_L_nReps100_nFolds4.csv'); 
WO_RF_Q = loadmodres('2016dec12_WO_RF_Q_nReps100_nFolds4.csv');

FM_LS_L = loadmodres('2016dec12_FM_LS_L_nReps100_nFolds4.csv'); 
FM_LS_Q = loadmodres('2016dec12_FM_LS_Q_nReps100_nFolds4.csv'); 

WO_LS_L = loadmodres('2016dec12_WO_LS_L_nReps100_nFolds4.csv'); 
WO_LS_Q = loadmodres('2016dec12_WO_LS_Q_nReps100_nFolds4.csv');

fid = fopen('yFM.csv'); pit = textscan(fid,'%s%s','delimiter',','); fclose(fid); yFM = str2double(pit{2}(2:end));
fid = fopen('yWO.csv'); pit = textscan(fid,'%s%s','delimiter',','); fclose(fid); yWO = str2double(pit{2}(2:end));

clearvars pit fid ans
%% Intialize variables to hold RMSE and R2
nreps = size(FM_LS_L,2);
RMSE = struct('FM_LS_L',zeros(1,nreps),'FM_LS_Q',zeros(1,nreps),'FM_RF_L',zeros(1,nreps),'FM_RF_Q',zeros(1,nreps),...
    'WO_LS_L',zeros(1,nreps),'WO_LS_Q',zeros(1,nreps),'WO_RF_L',zeros(1,nreps),'WO_RF_Q',zeros(1,nreps));
R2 = struct('FM_LS_L',zeros(1,nreps),'FM_LS_Q',zeros(1,nreps),'FM_RF_L',zeros(1,nreps),'FM_RF_Q',zeros(1,nreps),...
    'WO_LS_L',zeros(1,nreps),'WO_LS_Q',zeros(1,nreps),'WO_RF_L',zeros(1,nreps),'WO_RF_Q',zeros(1,nreps));
AR2 = struct('FM_LS_L',zeros(1,nreps),'FM_LS_Q',zeros(1,nreps),'FM_RF_L',zeros(1,nreps),'FM_RF_Q',zeros(1,nreps),...
    'WO_LS_L',zeros(1,nreps),'WO_LS_Q',zeros(1,nreps),'WO_RF_L',zeros(1,nreps),'WO_RF_Q',zeros(1,nreps));
SLP = struct('FM_LS_L',zeros(1,nreps),'FM_LS_Q',zeros(1,nreps),'FM_RF_L',zeros(1,nreps),'FM_RF_Q',zeros(1,nreps),...
    'WO_LS_L',zeros(1,nreps),'WO_LS_Q',zeros(1,nreps),'WO_RF_L',zeros(1,nreps),'WO_RF_Q',zeros(1,nreps));

%% Get all RMSE and Rsquared
for i = 1:nreps
    trash = fitlm(FM_LS_L(:,i),yFM); RMSE.FM_LS_L(i) = trash.RMSE; R2.FM_LS_L(i) = trash.Rsquared.Ordinary; AR2.FM_LS_L(i) = trash.Rsquared.Adjusted; SLP.FM_LS_L(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_LS_Q(:,i),yFM); RMSE.FM_LS_Q(i) = trash.RMSE; R2.FM_LS_Q(i) = trash.Rsquared.Ordinary; AR2.FM_LS_Q(i) = trash.Rsquared.Adjusted; SLP.FM_LS_Q(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_RF_L(:,i),yFM); RMSE.FM_RF_L(i) = trash.RMSE; R2.FM_RF_L(i) = trash.Rsquared.Ordinary; AR2.FM_RF_L(i) = trash.Rsquared.Adjusted; SLP.FM_RF_L(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(FM_RF_Q(:,i),yFM); RMSE.FM_RF_Q(i) = trash.RMSE; R2.FM_RF_Q(i) = trash.Rsquared.Ordinary; AR2.FM_RF_Q(i) = trash.Rsquared.Adjusted; SLP.FM_RF_Q(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_LS_L(:,i),yWO); RMSE.WO_LS_L(i) = trash.RMSE; R2.WO_LS_L(i) = trash.Rsquared.Ordinary; AR2.WO_LS_L(i) = trash.Rsquared.Adjusted; SLP.WO_LS_L(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_LS_Q(:,i),yWO); RMSE.WO_LS_Q(i) = trash.RMSE; R2.WO_LS_Q(i) = trash.Rsquared.Ordinary; AR2.WO_LS_Q(i) = trash.Rsquared.Adjusted; SLP.WO_LS_Q(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_RF_L(:,i),yWO); RMSE.WO_RF_L(i) = trash.RMSE; R2.WO_RF_L(i) = trash.Rsquared.Ordinary; AR2.WO_RF_L(i) = trash.Rsquared.Adjusted; SLP.WO_RF_L(i) = trash.Coefficients.Estimate(2);
    trash = fitlm(WO_RF_Q(:,i),yWO); RMSE.WO_RF_Q(i) = trash.RMSE; R2.WO_RF_Q(i) = trash.Rsquared.Ordinary; AR2.WO_RF_Q(i) = trash.Rsquared.Adjusted; SLP.WO_RF_Q(i) = trash.Coefficients.Estimate(2);
end

%% calculate mean and std for all model types

RMSE_M = zeros(1,8); RMSE_STD = zeros(1,8); R2_M = zeros(1,8); R2_STD = zeros(1,8);
AR2_M = zeros(1,8); AR2_STD = zeros(1,8); SLP_M = zeros(1,8); SLP_STD = zeros(1,8);

RMSE_M(1) = mean(RMSE.FM_LS_L); RMSE_M(2) = mean(RMSE.FM_LS_Q); RMSE_M(3) = mean(RMSE.FM_RF_L); RMSE_M(4) = mean(RMSE.FM_RF_Q);
RMSE_STD(1) = std(RMSE.FM_LS_L); RMSE_STD(2) = std(RMSE.FM_LS_Q); RMSE_STD(3) = std(RMSE.FM_RF_L); RMSE_STD(4) = std(RMSE.FM_RF_Q);

R2_M(1) = mean(R2.FM_LS_L); R2_M(2) = mean(R2.FM_LS_Q); R2_M(3) = mean(R2.FM_RF_L); R2_M(4) = mean(R2.FM_RF_Q);
R2_STD(1) = std(R2.FM_LS_L); R2_STD(2) = std(R2.FM_LS_Q); R2_STD(3) = std(R2.FM_RF_L); R2_STD(4) = std(R2.FM_RF_Q);

AR2_M(1) = mean(AR2.FM_LS_L); AR2_M(2) = mean(AR2.FM_LS_Q); AR2_M(3) = mean(AR2.FM_RF_L); AR2_M(4) = mean(AR2.FM_RF_Q);
AR2_STD(1) = std(AR2.FM_LS_L); AR2_STD(2) = std(AR2.FM_LS_Q); AR2_STD(3) = std(AR2.FM_RF_L); AR2_STD(4) = std(AR2.FM_RF_Q);

SLP_M(1) = mean(SLP.FM_LS_L); SLP_M(2) = mean(SLP.FM_LS_Q); SLP_M(3) = mean(SLP.FM_RF_L); SLP_M(4) = mean(SLP.FM_RF_Q);
SLP_STD(1) = std(SLP.FM_LS_L); SLP_STD(2) = std(SLP.FM_LS_Q); SLP_STD(3) = std(SLP.FM_RF_L); SLP_STD(4) = std(SLP.FM_RF_Q);

RMSE_M(5) = mean(RMSE.WO_LS_L); RMSE_M(6) = mean(RMSE.WO_LS_Q); RMSE_M(7) = mean(RMSE.WO_RF_L); RMSE_M(8) = mean(RMSE.WO_RF_Q);
RMSE_STD(5) = std(RMSE.WO_LS_L); RMSE_STD(6) = std(RMSE.WO_LS_Q); RMSE_STD(7) = std(RMSE.WO_RF_L); RMSE_STD(8) = std(RMSE.WO_RF_Q);

R2_M(5) = mean(R2.WO_LS_L); R2_M(6) = mean(R2.WO_LS_Q); R2_M(7) = mean(R2.WO_RF_L); R2_M(8) = mean(R2.WO_RF_Q);
R2_STD(5) = std(R2.WO_LS_L); R2_STD(6) = std(R2.WO_LS_Q); R2_STD(7) = std(R2.WO_RF_L); R2_STD(8) = std(R2.WO_RF_Q);

AR2_M(5) = mean(AR2.WO_LS_L); AR2_M(6) = mean(AR2.WO_LS_Q); AR2_M(7) = mean(AR2.WO_RF_L); AR2_M(8) = mean(AR2.WO_RF_Q);
AR2_STD(5) = std(AR2.WO_LS_L); AR2_STD(6) = std(AR2.WO_LS_Q); AR2_STD(7) = std(AR2.WO_RF_L); AR2_STD(8) = std(AR2.WO_RF_Q);

SLP_M(5) = mean(SLP.WO_LS_L); SLP_M(6) = mean(SLP.WO_LS_Q); SLP_M(7) = mean(SLP.WO_RF_L); SLP_M(8) = mean(SLP.WO_RF_Q);
SLP_STD(5) = std(SLP.WO_LS_L); SLP_STD(6) = std(SLP.WO_LS_Q); SLP_STD(7) = std(SLP.WO_RF_L); SLP_STD(8) = std(SLP.WO_RF_Q);

%% Let's get the bar graphs
% plot_comparison_bars(RMSE_M,RMSE_STD,'RMSE')
% plot_comparison_bars(R2_M,R2_STD,'Coefficient of Determination')
% plot_comparison_bars(AR2_M,AR2_STD,'Adjusted R2')
% plot_comparison_bars(SLP_M,SLP_STD,'Slope')

%% Finally let's plot both regular and adjusted R^2 in one plot
% col=hsv(11);
% col2 = col./2;
col= zeros(4,3);
col(1,:) = [1 0 0]; col(2,:) = [0 0 1]; col(3,:) = [1 0.7 0]; col(4,:) = [0 1 0];

figure
axes('Position',[0.05 0.02 0.95 0.95])
hold on

FullFMAR2Mat = [AR2.FM_LS_L;AR2.FM_LS_Q;AR2.FM_RF_L;AR2.FM_RF_Q];
FullWOAR2Mat = [AR2.WO_LS_L;AR2.WO_LS_Q;AR2.WO_RF_L;AR2.WO_RF_Q];

% inputs to this function: vector of things to plot, color of text inside
% the plot (make sure to check that against the variable Gr defined inside
% the function), xlimits for x-plotting, color of dots plotted, and scale
% factor for the x-scatter.
scf = 50; % scale factor for x-scatter
plot1DDistributionV2(FullFMAR2Mat(1,:),'k',[0.9,1.1],col(1,:),scf) 
plot1DDistributionV2(FullFMAR2Mat(2,:),'k',[1.9,2.1],col(2,:),scf)
plot1DDistributionV2(FullFMAR2Mat(3,:),'k',[2.9,3.1],col(3,:),scf)
plot1DDistributionV2(FullFMAR2Mat(4,:),'k',[3.9,4.1],col(4,:),scf)

plot1DDistributionV2(FullWOAR2Mat(1,:),'k',[5.9,6.1],col(1,:),scf)
plot1DDistributionV2(FullWOAR2Mat(2,:),'k',[6.9,7.1],col(2,:),scf)
plot1DDistributionV2(FullWOAR2Mat(3,:),'k',[7.9,8.1],col(3,:),scf)
plot1DDistributionV2(FullWOAR2Mat(4,:),'k',[8.9,9.1],col(4,:),scf)

ax = gca;
% ax.XTick = [2.5 7.5];
ax.XTick = [];
% ax.TickLength = [0 0];
% ax.XTickLabel = {'Fugl-Meyer','Wolf Motor Function'};
set(ax,'FontSize',18)
ylabel('Adjusted Coefficient of Determination R^2')
% ylabel('R^2 (light)/ Adjusted R^2 (dark)')
ylim([0 1.1])
xlim([0 10])
% text(1,0.95,'First order lasso','Rotation',90,'FontSize',21,'HorizontalAlignment','left')
% text(2,0.95,'Second order lasso','Rotation',90,'FontSize',21,'HorizontalAlignment','left')
% text(3,0.95,'First order random forests','Rotation',90,'FontSize',21,'HorizontalAlignment','left')
% text(4,0.95,'Second order random forests','Rotation',90,'FontSize',21,'HorizontalAlignment','left')
% 
% text(6,0.95,'First order lasso','Rotation',90,'FontSize',21,'HorizontalAlignment','left')
% text(7,0.95,'Second order lasso','Rotation',90,'FontSize',21,'HorizontalAlignment','left')
% text(8,0.95,'First order random forests','Rotation',90,'FontSize',21,'HorizontalAlignment','left')
% text(9,0.95,'Second order random forests','Rotation',90,'FontSize',21,'HorizontalAlignment','left')


