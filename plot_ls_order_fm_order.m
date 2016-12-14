
clear all
clc
%% This script will do a simple plot of percentage of times a feature was selected by lasso, and also the 

[FM_order_data,~,~] = xlsread('LS_FM_order.csv');
[WO_order_data,~,~] = xlsread('LS_WO_order.csv');

[FM_VIs_data,Rownames,~] = xlsread('RF_VIs_FM_Raw_v8.csv');
[WO_VIs_data,~,~] = xlsread('RF_VIs_WO_Raw_v8.csv');
Rownames(1,:) = [];
Rownames(:,2:end) = [];

%% figure out the order for VIs, also normalize VIs

% [~,FM_VI_order] = sort(mean(FM_VIs_data,2),'descend');
% [~,WO_VI_order] = sort(mean(WO_VIs_data,2),'descend');

% FM_VI_norm = (FM_VIs_data)./(max(max(FM_VIs_data)));
% WO_VI_norm = (WO_VIs_data)./(max(max(WO_VIs_data)));

[~,FM_LS_order] = sort(FM_order_data,'descend');
[~,WO_LS_order] = sort(WO_order_data,'descend');

%% let's plot

figure
clf
xx = linspace(15,1,15);
subplot(1,2,1);hold on;
% for i = 1:15
%     plot(FM_VI_norm(FM_VI_order(i),:),(15-i+1)*ones(400,1),'b.')
% end
plot(FM_order_data(FM_LS_order(1:15)),xx,'bd','MarkerSize',7,'MarkerFaceColor','b')
xlim([0 1.1])
ylim([0.5 15.5])
for i = 1:15
    text(0.2,15-i+1,Rownames(FM_LS_order(i)));
end

subplot(1,2,2);hold on;
plot(WO_order_data(WO_LS_order(1:15)),xx,'bd','MarkerSize',7,'MarkerFaceColor','b')
xlim([0 1.1])
ylim([0.5 15.5])
for i = 1:15
    text(0.2,15-i+1,Rownames(WO_LS_order(i)));
end

%% plot top features
[DF1, txt_DF1, ~] = xlsread('df1.csv');
yFM = xlsread('yFM.csv');
yWO = xlsread('yWO.csv');
Names_df1 = txt_DF1(1,:);

DF1(:,1) = [];
yFM(:,1) = [];
yWO(:,1) = [];

figure
clf
% from lasso used to predict
% top5FM = FM_LS_order(1:5);
% top5WO = WO_LS_order(1:5);
% from RF used to affect
top5FM = [5,42,39,9,18];
top5WO = [53,54,24,12,52];

for i = 1:5
    subplot(2,5,i); X = [ones(26,1),DF1(:,top5FM(i))]; b = X\yFM; CalcFM = X*b;
    scatter(DF1(:,top5FM(i)),yFM,'MarkerFaceColor','b','MarkerEdgeColor','b');
    hold on
    plot(DF1(:,top5FM(i)),CalcFM,'r-','LineWidth',1)
    xlabel(Names_df1(top5FM(i)))
    if i == 1
        ylabel('Change in Fugl-Meyer')
    else
        ylabel('')
    end
    if max(DF1(:,top5FM(i))) == 1
        xlim([-0.2 1.2])
    end
    
    subplot(2,5,i+5); X = [ones(26,1),DF1(:,top5WO(i))]; b = X\yWO; CalcWO = X*b;
    scatter(DF1(:,top5WO(i)),yWO,'MarkerFaceColor','b','MarkerEdgeColor','b');
    hold on
    plot(DF1(:,top5WO(i)),CalcWO,'r-','LineWidth',1)
    xlabel(Names_df1(top5WO(i)))
    if i == 1
        ylabel('Change in Wolf Motor Function')
    else
        ylabel('')
    end
    if max(DF1(:,top5WO(i))) == 1
        xlim([-0.2 1.2])
    end
    
end
