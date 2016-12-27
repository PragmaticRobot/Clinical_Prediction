% This script imports R results for bottom up models (lasso_into_lasso and RF_into_lasso)

clear all
clc

%% load the needed data and clean up (first column is row number, should be removed

% for 26 choose 1 (LOO-CV)
WO_RF_L = xlsread('2016dec12_WO_RF_L_R2_meta_26C1.csv');
WO_RF_Q = xlsread('2016dec12_WO_RF_Q_R2_meta_26C1.csv');
WO_LS_L = xlsread('2016dec12_WO_LS_Q_R2_meta_26C1.csv');
WO_LS_Q = xlsread('2016dec12_WO_LS_L_R2_meta_26C1.csv');

% for 26 choose 2 (leave 2 out)
% WO_RF_L = xlsread('WO_RF_L_26C2.csv');
% WO_RF_Q = xlsread('WO_RF_Q_26C2.csv');
% WO_LS_L = xlsread('WO_LS_L_26C2.csv');
% WO_LS_Q = xlsread('WO_LS_Q_26C2.csv');

WO_RF_L(:,1) = []; WO_RF_Q(:,1) = []; WO_LS_L(:,1) = []; WO_LS_Q(:,1) = [];

%% Define vertices for patch objects:
V1 = 1:20;
V2 = linspace(20,1,20);
Vx = horzcat(V1,V2);

V_WO_LS_L = mean(WO_LS_L,1) - std(WO_LS_L,1);
V_WO_LS_Q = mean(WO_LS_Q,1) - std(WO_LS_Q,1);
V_WO_RF_L = mean(WO_RF_L,1) - std(WO_RF_L,1);
V_WO_RF_Q = mean(WO_RF_Q,1) - std(WO_RF_Q,1);

V_WO_LS_L2 = mean(WO_LS_L,1) + std(WO_LS_L,1);
V_WO_LS_Q2 = mean(WO_LS_Q,1) + std(WO_LS_Q,1);
V_WO_RF_L2 = mean(WO_RF_L,1) + std(WO_RF_L,1);
V_WO_RF_Q2 = mean(WO_RF_Q,1) + std(WO_RF_Q,1);

VL1 = [V_WO_LS_L,V_WO_LS_L2(end:-1:1)];
VL2 = [V_WO_LS_Q,V_WO_LS_Q2(end:-1:1)];
VL3 = [V_WO_RF_L,V_WO_RF_L2(end:-1:1)];
VL4 = [V_WO_RF_Q,V_WO_RF_Q2(end:-1:1)];

%% Now the plotting

colors(1,:) = [1 0 0]; colors(2,:) = [0 0 1]; colors(3,:) = [1 0.7 0]; colors(4,:) = [0 1 0];

% colors(1,:) = [0.933 0.172 0.172]; colors(2,:) = [0 0.545 0.545]; colors(3,:) = [1 0.5 0]; colors(4,:) = [110, 152, 135]/255;

xx2 = 1:20;
figure
clf
hold on

% xx = repmat(1:20,26,1);
xx = repmat(1:20,325,1); % for leave 2 out
% for i = 1:20;
%     scatter(xx(:,i),WO_LS_L(:,i),1,colors(1,:))
%     scatter(xx(:,i),WO_LS_Q(:,i),1,colors(2,:))
%     scatter(xx(:,i),WO_RF_L(:,i),1,colors(3,:))
%     scatter(xx(:,i),WO_RF_Q(:,i),1,colors(4,:))
% end

q1 = plot(xx2,mean(WO_LS_L,1),'Color',colors(1,:),'LineWidth',1.5);
q2 = plot(xx2,mean(WO_LS_Q,1),'Color',colors(2,:),'LineWidth',1.5);
q3 = plot(xx2,mean(WO_RF_L,1),'Color',colors(3,:),'LineWidth',1.5);
q4 = plot(xx2,mean(WO_RF_Q,1),'Color',colors(4,:),'LineWidth',1.5);

patch(Vx,VL1,colors(1,:),'FaceAlpha',0.3,'EdgeColor','none');
patch(Vx,VL2,colors(2,:),'FaceAlpha',0.3,'EdgeColor','none');
patch(Vx,VL3,colors(3,:),'FaceAlpha',0.3,'EdgeColor','none');
patch(Vx,VL4,colors(4,:),'FaceAlpha',0.3,'EdgeColor','none');

ax = gca;
ax.XTick = 1:1:20;
xlabel('Number of base features')
ylim([0 1.1])
legend([q1;q2;q3;q4],{'1st order lasso';'2nd order lasso';'1st order RF';'2nd order RF'})