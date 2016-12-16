clear all
clc
colormap winter
% load the lasso counts results

[FM_LS, txt_FM_LS,~] = xlsread('2016dec12_FM_lasso_counts.csv');
[WO_LS, ~,~] = xlsread('2016dec12_WO_lasso_counts.csv');
featNames = txt_FM_LS(2:end,1);
colnames = txt_FM_LS(1,2:end);

FM_counts = table(FM_LS(:,1),FM_LS(:,2),FM_LS(:,3), 'VariableNames',colnames,'RowNames',featNames);
WO_counts = table(WO_LS(:,1),WO_LS(:,2),WO_LS(:,3), 'VariableNames',colnames,'RowNames',featNames);

% clearvars -except FM_counts WO_counts

bar(FM_LS)
legend('4 fold good','LOO good','13 fold good')
ax = gca;
ax.XTick = 1:1:52;
ax.XTickLabel = featNames;
ax.XTickLabelRotation = 90;
title('lasso features picked by different cross-validations')

%% let plot relationship between features and outcome for "winners"
[DF1, txt_DF1, ~] = xlsread('df1.csv');
[DF3, txt_DF3, ~] = xlsread('df3.csv');
yFM = xlsread('yFM.csv');
yWO = xlsread('yWO.csv');
Names_df1 = txt_DF1(1,:);
Names_df3 = txt_DF3(1,:);
DF1(:,1) = [];
DF3(:,1) = [];
yFM(:,1) = [];
yWO(:,1) = [];

%% plot the top 5 for each
% from lasso used to predict
% top5FM = [39,46,33,27,40];
% top5WO = [53,45,40,38,19];
% from RF used to affect
top5FM = [5,33,42,39,18];
top5WO = [53,54,52,12,24];

figure
for i = 1:5
    subplot(2,5,i); X = [ones(26,1),DF1(:,top5FM(i))]; b = X\yFM; CalcFM = X*b;
    scatter(DF1(:,top5FM(i)),yFM);
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
    scatter(DF1(:,top5WO(i)),yWO);
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