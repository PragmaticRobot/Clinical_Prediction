% This script takes the results of the R models (LASSO, RF) with ranks for
% both Fugl-Meyer (FM) and Wolf (WO).
clear all
clc

%% data prep
[FM_VI_data,RownamesFM,~] = xlsread('2016dec12_RF_Best_FM_VI.csv');
[WO_VI_data,RownamesWO,~] = xlsread('2016dec12_RF_Best_WO_VI.csv');

[~,FM_RF_order] = sort(mean(FM_VI_data,2),'descend');
[~,WO_RF_order] = sort(mean(WO_VI_data,2),'descend');

[DF1, txt_DF1, ~] = xlsread('2016dec12_df1.csv'); DF1(:,1) = [];
yFM = xlsread('yFM.csv'); yFM(:,1) = [];
yWO = xlsread('yWO.csv'); yWO(:,1) = [];
Names_df1 = txt_DF1(1,:);
top5FM = FM_RF_order(1:5);
top5WO = WO_RF_order(1:5);

% import lasso results
[FM_LS_counts,~,~] = xlsread('2016dec12_LS_FM_order.csv');
[WO_LS_counts,~,~] = xlsread('2016dec12_LS_WO_order.csv');

%% Cleanup headings for name fields
RownamesFM(1,:) = [];
RownamesFM(:,2:end) = [];

RownamesWO(1,:) = [];
RownamesWO(:,2:end) = [];

%% convert data type so it's more usable

FM_RF = FM_VI_data;
WO_RF = WO_VI_data;

%% plotting Random Forests Fugl-Meyer

% box limits
Y = prctile(FM_RF,[25 50 75],2); % calculate 25% 50% 75% for each feature
xmax = ceil(max(max(FM_RF)));
xmin = floor(min(min(FM_RF)));
Y = Y/xmax;
nrow = size(FM_RF,1);

fig = figure;
fig.PaperType = 'usletter';
fig.PaperOrientation = 'portrait';
clf
% subplot(1,2,1);
axes('Position',[0.04 .05 .44 .94])
hold on
ylim([0 nrow+0.5])
% text(0.84,47.6,RownamesFM(1),'HorizontalAlignment','right');

for i = 1:nrow
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;
    
    % add the dot from lasso, fraction of cross-validations where
    % this feature was used
    plot(FM_LS_counts(FM_RF_order(i)),ymid,'bd','MarkerSize',7,'MarkerFaceColor',[92,171,181]/255,'MarkerEdgeColor','none');
    
    % plot RF ranks
    for j = 1:100
        yy = ymid + ((rand-0.5)*.7); % add noise to y-value
        plot(FM_RF(FM_RF_order(i),j)/xmax,yy,'.','Color',0.2*[1 1 1],'MarkerSize',4);
    end
    % now plot the box for each value
    line([Y(FM_RF_order(i),1) Y(FM_RF_order(i),3)],[ymid-0.25 ymid-0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(FM_RF_order(i),1) Y(FM_RF_order(i),1)],[ymid-0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(FM_RF_order(i),3) Y(FM_RF_order(i),3)],[ymid-0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(FM_RF_order(i),3) Y(FM_RF_order(i),1)],[ymid+0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(FM_RF_order(i),2) Y(FM_RF_order(i),2)],[ymid-0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1); % Median

    % text with the feature name
    if i == 1
        text(Y(FM_RF_order(i),2)-0.03, ymid, RownamesFM(FM_RF_order(i)),'HorizontalAlignment','right');
    else
        text(Y(FM_RF_order(i),2)+0.03, ymid, RownamesFM(FM_RF_order(i)),'HorizontalAlignment','left');
    end

end
box off
set(gca,'color','none','YColor','none')%,'Ticklength', [0 0])
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1])
xlim([xmin/xmax 1])
xlabel('Normalized Feature Importance')
% draw the extra axis for lasso
line([0 1],[51 51],'Color',[92,171,181]/255); 
line([1 1],[51 51.3],'Color',[92,171,181]/255);        text(1,51.5,'1','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.8 0.8],[51 51.3],'Color',[92,171,181]/255);    text(0.8,51.5,'0.8','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.6 0.6],[51 51.3],'Color',[92,171,181]/255);    text(0.6,51.5,'0.6','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.4 0.4],[51 51.3],'Color',[92,171,181]/255);    text(0.4,51.5,'0.4','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.2 0.2],[51 51.3],'Color',[92,171,181]/255);    text(0.2,51.5,'0.2','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0 0],[51 51.3],'Color',[92,171,181]/255);        text(0,51.5,'0','Color',[92,171,181]/255,'HorizontalAlignment','center');

% scale yFM for the smaller plots, divide by max yFM
maxFM = max(yFM);
yFM_sc = yFM./maxFM; % this scales it between -1 and 1
yFM_sc = yFM_sc./2; % this scales it between -0.5 and 0.5

% plot individual small plots
axes('Position',[0.02 .05 .01 0.94])
hold on
ylim([0 nrow+0.5])
xlim([-0.5 0.8])
for i = 1:nrow
    ymid = (nrow-i)+0.5;
    % plot individual top features
    this_feat = DF1(:,FM_RF_order(i)); % the feature to be plotted at position i (counting from the top)
    this_feat_sc = this_feat./max(this_feat); % scale the feature between -1 and 1
    X = [ones(26,1),this_feat]; b = X\yFM; CalcFM = X*b;
    scatter(this_feat_sc-0.25,yFM_sc+ymid,1,'.','MarkerFaceColor','b','MarkerEdgeColor','b');
    hold on
    CalcFM_sc = (CalcFM./maxFM)./2;
    plot(this_feat_sc-0.25,CalcFM_sc+ymid,'Color','m','LineWidth',1)
end
axis off
box off
set(gca,'color','none','YColor','none','Ticklength', [0 0])
ax = gca; ax.XTick = []; ax.YTick = [];

%% plotting Random Forests Wolf

% box limits
Y = prctile(WO_RF,[25 50 75],2); % calculate 25% 50% 75% for each feature
xmax = ceil(max(max(WO_RF)));
% xmin = floor(min(min(WO_RF(1:15,:))));
xmin = floor(min(min(WO_RF)));
Y = Y/xmax;

axes('Position',[.54 .05 .43 .94])
hold on
ylim([0 nrow+0.5])

for i = 1:nrow
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;
    % add the dot from lasso, fraction of cross-validations where
    % this feature was used
    plot(WO_LS_counts(WO_RF_order(i)),ymid,'bd','MarkerSize',7,'MarkerFaceColor',[92,171,181]/255,'MarkerEdgeColor','none');
    % now plot RF
    for j = 1:100
        yy = ymid + ((rand-0.5)*.7); % add noise to y-value
        plot(WO_RF(WO_RF_order(i),j)/xmax,yy,'.','Color',0.2*[1 1 1],'MarkerSize',4);
    end
%     text(xmin/xmax*8,ymid,rnam4{i},'HorizontalAlignment','left');
    % now plot the box for each value
    line([Y(WO_RF_order(i),1) Y(WO_RF_order(i),3)],[ymid-0.25 ymid-0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(WO_RF_order(i),1) Y(WO_RF_order(i),1)],[ymid-0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(WO_RF_order(i),3) Y(WO_RF_order(i),3)],[ymid-0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(WO_RF_order(i),3) Y(WO_RF_order(i),1)],[ymid+0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1);
    line([Y(WO_RF_order(i),2) Y(WO_RF_order(i),2)],[ymid-0.25 ymid+0.25],'Color',0.8*[1 0.7 0],'LineWidth',1); % Median
    
    if i == 1
        text(Y(WO_RF_order(i),2)-0.03, ymid, RownamesWO(WO_RF_order(i)),'HorizontalAlignment','right');
    else
        text(Y(WO_RF_order(i),2)+0.03, ymid, RownamesWO(WO_RF_order(i)),'HorizontalAlignment','left');
    end
end
box off
set(gca,'color','none','YColor','none')%,'Ticklength', [0 0])
set(gca,'XTick',[0 0.2 0.4 0.6 0.8 1])
xlim([xmin/xmax 1])
xlabel('Normalized Feature Importance')

% draw the extra axis for lasso
line([0 1],[51 51],'Color',[92,171,181]/255); 
line([1 1],[51 51.3],'Color',[92,171,181]/255);        text(1,51.5,'1','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.8 0.8],[51 51.3],'Color',[92,171,181]/255);    text(0.8,51.5,'0.8','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.6 0.6],[51 51.3],'Color',[92,171,181]/255);    text(0.6,51.5,'0.6','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.4 0.4],[51 51.3],'Color',[92,171,181]/255);    text(0.4,51.5,'0.4','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0.2 0.2],[51 51.3],'Color',[92,171,181]/255);    text(0.2,51.5,'0.2','Color',[92,171,181]/255,'HorizontalAlignment','center');
line([0 0],[51 51.3],'Color',[92,171,181]/255);        text(0,51.5,'0','Color',[92,171,181]/255,'HorizontalAlignment','center');

% scale yFM for the smaller plots, divide by max yFM
maxWO = max(yWO);
yWO_sc = yWO./maxWO; % this scales it between -1 and 1
yWO_sc = yWO_sc./2; % this scales it between -0.5 and 0.5

% now plot each feature
axes('Position',[0.52 .05 .01 0.94])
hold on
ylim([0 nrow+0.5])
xlim([-0.5 0.8])
for i = 1:nrow
    ymid = (nrow-i)+0.5;
    % plot individual top features
    this_feat = DF1(:,WO_RF_order(i)); % the feature to be plotted at position i (counting from the top)
    this_feat_sc = this_feat./max(this_feat); % scale the feature between -1 and 1
    X = [ones(26,1),this_feat]; b = X\yWO; CalcWO = X*b;
    scatter(this_feat_sc-0.25,yWO_sc+ymid,1,'.','MarkerFaceColor','b','MarkerEdgeColor','b');
    hold on
    CalcWO_sc = (CalcWO./maxWO)./2;
    plot(this_feat_sc-0.25,CalcWO_sc+ymid,'Color','m','LineWidth',1)
end
axis off
box off
set(gca,'color','none','YColor','none','Ticklength', [0 0])
ax = gca; ax.XTick = []; ax.YTick = [];
