clear all
clc
%% Load all the files
% the models are in order: Lin, NoExtrap, Quad
% 1-2: FM LS
% 4-5: FM RF
% 7-8: WO LS
% 10-11: WO RF
FM_LS_L = loadmodres('mod1_v11.csv'); FM_LS_Q = loadmodres('mod2_v11.csv'); FM_RF_L = loadmodres('mod3_v8.csv'); FM_RF_Q = loadmodres('mod4_v8.csv');

WO_LS_L = loadmodres('mod5_v11.csv'); WO_LS_Q = loadmodres('mod6_v11.csv'); WO_RF_L = loadmodres('mod7_v8.csv'); WO_RF_Q = loadmodres('mod8_v8.csv');

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
plot_comparison_bars(RMSE_M,RMSE_STD,'RMSE')
plot_comparison_bars(R2_M,R2_STD,'Coefficient of Determination')
plot_comparison_bars(AR2_M,AR2_STD,'Adjusted R2')
plot_comparison_bars(SLP_M,SLP_STD,'Slope')

%% Finally let's plot both regular and adjusted R^2 in one plot
% col=hsv(11);
% col2 = col./2;
col= zeros(4,3);
col(1,:) = [1 0 0]; col(2,:) = [0 0 1]; col(3,:) = [1 0.7 0]; col(4,:) = [0 1 0];
col2 = col./2;

figure
hold on

xxx = 1:1:9;
FullFMR2Mat = [R2.FM_LS_L;R2.FM_LS_Q;R2.FM_RF_L;R2.FM_RF_Q];
FullWOR2Mat = [R2.WO_LS_L;R2.WO_LS_Q;R2.WO_RF_L;R2.WO_RF_Q];
FullFMAR2Mat = [AR2.FM_LS_L;AR2.FM_LS_Q;AR2.FM_RF_L;AR2.FM_RF_Q];
FullWOAR2Mat = [AR2.WO_LS_L;AR2.WO_LS_Q;AR2.WO_RF_L;AR2.WO_RF_Q];
X = (rand(100,1)-0.5)/10;

plot(1.0*ones(1,100)+X',FullFMR2Mat(1,:),'.','color',col(1,:),'MarkerSize',7);
plot((1+0.2).*ones(1,100)+X',FullFMAR2Mat(1,:),'.','color',col2(1,:),'MarkerSize',7);
plot(2.0*ones(1,100)+X',FullFMR2Mat(2,:),'.','color',col(2,:),'MarkerSize',7);
plot((2+0.2).*ones(1,100)+X',FullFMAR2Mat(2,:),'.','color',col2(2,:),'MarkerSize',7);
plot(3.0*ones(1,100)+X',FullFMR2Mat(3,:),'.','color',col(3,:),'MarkerSize',7);
plot((3+0.2).*ones(1,100)+X',FullFMAR2Mat(3,:),'.','color',col2(3,:),'MarkerSize',7);
plot(4.0*ones(1,100)+X',FullFMR2Mat(4,:),'.','color',col(4,:),'MarkerSize',7);
plot((4+0.2).*ones(1,100)+X',FullFMAR2Mat(4,:),'.','color',col2(4,:),'MarkerSize',7);

plot(6.0*ones(1,100)+X',FullWOR2Mat(1,:),'.','color',col(1,:),'MarkerSize',7);
plot((6+0.2).*ones(1,100)+X',FullWOAR2Mat(1,:),'.','color',col2(1,:),'MarkerSize',7);
plot(7.0*ones(1,100)+X',FullWOR2Mat(2,:),'.','color',col(2,:),'MarkerSize',7);
plot((7+0.2).*ones(1,100)+X',FullWOAR2Mat(2,:),'.','color',col2(2,:),'MarkerSize',7);
plot(8.0*ones(1,100)+X',FullWOR2Mat(3,:),'.','color',col(3,:),'MarkerSize',7);
plot((8+0.2).*ones(1,100)+X',FullWOAR2Mat(3,:),'.','color',col2(3,:),'MarkerSize',7);
plot(9.0*ones(1,100)+X',FullWOR2Mat(4,:),'.','color',col(4,:),'MarkerSize',7);
plot((9+0.2).*ones(1,100)+X',FullWOAR2Mat(4,:),'.','color',col2(4,:),'MarkerSize',7);

for i = 1:8
    if i <= 4
        plot(xxx(i),R2_M(i),'d','MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2);
        plot(xxx(i)+0.2,AR2_M(i),'d','MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2);
        plot( [xxx(i) xxx(i)], R2_M(i)+R2_STD(i)*[1 -1], 'color','k',  'lineWidth',2);
        plot( [xxx(i)+0.2 xxx(i)+0.2], AR2_M(i)+AR2_STD(i)*[1 -1], 'color','k',  'lineWidth',2);
       
    else
        j = i+1;
        plot(xxx(j),R2_M(i),'d','MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2);
        plot(xxx(j)+0.2,AR2_M(i),'d','MarkerEdgeColor','k','MarkerSize',12,'LineWidth',2);
        plot( [xxx(j) xxx(j)], R2_M(i)+R2_STD(i)*[1 -1], 'color','k',  'lineWidth',2);
        plot( [xxx(j)+0.2 xxx(j)+0.2], AR2_M(i)+AR2_STD(i)*[1 -1], 'color','k',  'lineWidth',2);
    end
end


ax = gca;
ax.XTick = [2.5 7.5];
ax.XTickLabel = {'Fugl-Meyer','Wolf Motor Function'};
set(ax,'FontSize',18)
% ylabel('R^2 (light)/ Adjusted R^2 (dark)')
ylim([0 1.2])
xlim([0 10])
text(1,0.81,'First order lasso','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(2,0.84,'Second order lasso','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(3,0.15,'First order random forests','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(4,0.15,'Second order random forests','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(0.975,0.02,'R^2','color','k','Rotation',90,'FontSize',17)
text(1.175,0.02,'Adjusted R^2','color','k','Rotation',90,'FontSize',17)
text(1.975,0.02,'R^2','color','k','Rotation',90,'FontSize',17)
text(2.175,0.02,'Adjusted R^2','color','k','Rotation',90,'FontSize',17)


text(6,0.95,'First order lasso','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(7,0.9,'Second order lasso','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(8,0.15,'First order random forests','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(9,0.15,'Second order random forests','Rotation',60,'FontSize',17,'HorizontalAlignment','left')
text(5.975,0.02,'R^2','color','k','Rotation',90,'FontSize',17)
text(6.175,0.02,'Adjusted R^2','color','k','Rotation',90,'FontSize',17)
text(6.975,0.02,'R^2','color','k','Rotation',90,'FontSize',17)
text(7.175,0.02,'Adjusted R^2','color','k','Rotation',90,'FontSize',17)


