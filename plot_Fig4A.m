% This script takes the results of the R models (LASSO, RF) with ranks for
% both Fugl-Meyer (FM) and Wolf (WO).
clear all
clc

%% data prep
[FM_VI_data,RownamesFM,~] = xlsread('2016dec12_FM_RF_VI_sorted.csv');
[WO_VI_data,RownamesWO,~] = xlsread('2016dec12_WO_RF_VI_sorted.csv');

[FM_RF_order,~,~] = xlsread('2016dec12_FM_RF_Order.csv');
[WO_RF_order,~,~] = xlsread('2016dec12_WO_RF_Order.csv');

[DF1, txt_DF1, ~] = xlsread('2016dec12_df1.csv'); DF1(:,1) = [];
yFM = xlsread('yFM.csv'); yFM(:,1) = [];
yWO = xlsread('yWO.csv'); yWO(:,1) = [];
Names_df1 = txt_DF1(1,:);
top5FM = FM_RF_order(1:5,2);
top5WO = WO_RF_order(1:5,2);

%% Cleanup headings for name fields
RownamesFM(1,:) = [];
RownamesFM(:,2:end) = [];

RownamesWO(1,:) = [];
RownamesWO(:,2:end) = [];

%% convert data type so it's more usable

FM_RF = FM_VI_data;
WO_RF = WO_VI_data;

%% colormaps
% Frequency of selection will be conveyed through patches behind the plot,
% with the darker patches correlating with more frequent selection. so
% first we need to create a colormap:

% start values are the base values for colors here
% Red for Fugl-Meyer (matching R plots)
% start 238 59 59
% end 251 208 208
map = zeros(100,3);
map(:,1) = linspace(251/255,238/255,100);
map(:,2) = linspace(208/255,150/255,100);
map(:,3) = linspace(208/255,150/255,100);
mapComp = ones(100,3)-map;
% Green for Wolf (to match R plots)
% start 85 107 47
% end 232 239 220
map2 = zeros(100,3);
map2(:,1) = linspace(232/255,176/255,100);
map2(:,2) = linspace(222/255,200/255,100);
map2(:,3) = linspace(220/255,97/255,100);
map2Comp = ones(100,3)-map2;

col = 0.7*[1 0 0];
col2 = 0.7*[0 0 1];

%% plotting Random Forests Fugl-Meyer

% box limits
Y = prctile(FM_RF,[25 50 75],2); % calculate 25% 50% 75% for each feature
xmax = ceil(max(max(FM_RF)));
xmin = floor(min(min(FM_RF)));
Y = Y/xmax;

figure
clf
subplot(1,2,1);
set(gca,'Position',[0.05 .1 .45 .9])
hold on

% box around each feature
line([0.85 0.94],[46.5 46.5],'Color',0.2*[1 1 1]);
line([0.94 0.94],[46.5 48.8],'Color',0.2*[1 1 1]);
line([0.94 0.85],[48.8 48.8],'Color',0.2*[1 1 1]);
line([0.85 0.85],[48.8 46.5],'Color',0.2*[1 1 1]);

line([0.85 0.94],[46.3 46.3],'Color',0.2*[1 1 1]);
line([0.94 0.94],[46.3 44.0],'Color',0.2*[1 1 1]);
line([0.94 0.85],[44.0 44.0],'Color',0.2*[1 1 1]);
line([0.85 0.85],[44.0 46.3],'Color',0.2*[1 1 1]);

line([0.85 0.94],[43.8 43.8],'Color',0.2*[1 1 1]);
line([0.94 0.94],[43.8 41.5],'Color',0.2*[1 1 1]);
line([0.94 0.85],[41.5 41.5],'Color',0.2*[1 1 1]);
line([0.85 0.85],[41.5 43.8],'Color',0.2*[1 1 1]);

line([0.85 0.94],[41.3 41.3],'Color',0.2*[1 1 1]);
line([0.94 0.94],[41.3 39.0],'Color',0.2*[1 1 1]);
line([0.94 0.85],[39.0 39.0],'Color',0.2*[1 1 1]);
line([0.85 0.85],[39.0 41.3],'Color',0.2*[1 1 1]);

line([0.85 0.94],[38.8 38.8],'Color',0.2*[1 1 1]);
line([0.94 0.94],[38.8 36.5],'Color',0.2*[1 1 1]);
line([0.94 0.85],[36.5 36.5],'Color',0.2*[1 1 1]);
line([0.85 0.85],[36.5 38.8],'Color',0.2*[1 1 1]);

% surrounding box
line([0.7 0.95],[36.3 36.3],'Color',0.2*[1 1 1]);
line([0.95 0.95],[36.3 49.0],'Color',0.2*[1 1 1]);
line([0.95 0.7],[49.0 49.0],'Color',0.2*[1 1 1]);
line([0.7 0.7],[49.0 36.3],'Color',0.2*[1 1 1]);

text(0.84,47.6,RownamesFM(1),'HorizontalAlignment','right');
text(0.84,45.1,RownamesFM(2),'HorizontalAlignment','right');
text(0.84,42.6,RownamesFM(3),'HorizontalAlignment','right');
text(0.84,40.1,RownamesFM(4),'HorizontalAlignment','right');
text(0.84,37.6,RownamesFM(5),'HorizontalAlignment','right');

for i = 1:51
    nrow = size(FM_RF,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;
    for j = 1:100
        yy = ymid + ((rand-0.5)*.7); % add noise to y-value
        plot(FM_RF(i,j)/xmax,yy,'.','Color',[0 139/255 139/255],'MarkerSize',5);
    end
%     text(xmin/xmax*2.1,ymid,rnam2{i},'HorizontalAlignment','left');
    % now plot the box for each value
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',[205/255 102/255 0]);
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0]);
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0]);
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',[205/255 102/255 0]);
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0]); % Median
    if i >= 9
        text(Y(i,2)+0.03, ymid, RownamesFM(i),'HorizontalAlignment','left');
    else
        text(Y(i,2)-0.03, ymid, RownamesFM(i),'HorizontalAlignment','right');
    end
end
box off
set(gca,'color','none','YColor','none','Ticklength', [0 0])
xlim([xmin/xmax 1])
xlabel('Normalized Feature Importance')
% title('Variable Importance for Change in Fugl-Meyer')

% plot individual top features
axes('Position',[0.395 0.7 0.062 0.11])
X = [ones(26,1),DF1(:,top5FM(1))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(1)),yFM,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(1)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.413 0.56 0.065 0.11])
X = [ones(26,1),DF1(:,top5FM(2))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(2)),yFM,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(2)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.396 0.42 0.0654 0.11])
X = [ones(26,1),DF1(:,top5FM(3))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(3)),yFM,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(3)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.401 0.28 0.0654 0.11])
X = [ones(26,1),DF1(:,top5FM(4))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(4)),yFM,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(4)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.415 0.14 0.045 0.11])
X = [ones(26,1),DF1(:,top5FM(5))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(5)),yFM,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(5)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

%% plotting Random Forests Wolf

% box limits
Y = prctile(WO_RF,[25 50 75],2); % calculate 25% 50% 75% for each feature
xmax = ceil(max(max(WO_RF)));
% xmin = floor(min(min(WO_RF(1:15,:))));
xmin = floor(min(min(WO_RF)));
Y = Y/xmax;
% figure(2)
% clf
subplot(1,2,2);
set(gca,'Position',[.55 .1 .4 .9])
hold on
% box around each feature
line([0.80 0.94],[46.5 46.5],'Color',0.2*[1 1 1]);
line([0.94 0.94],[46.5 48.8],'Color',0.2*[1 1 1]);
line([0.94 0.80],[48.8 48.8],'Color',0.2*[1 1 1]);
line([0.80 0.80],[48.8 46.5],'Color',0.2*[1 1 1]);

line([0.80 0.94],[46.3 46.3],'Color',0.2*[1 1 1]);
line([0.94 0.94],[46.3 44.0],'Color',0.2*[1 1 1]);
line([0.94 0.80],[44.0 44.0],'Color',0.2*[1 1 1]);
line([0.80 0.80],[44.0 46.3],'Color',0.2*[1 1 1]);

line([0.80 0.94],[43.8 43.8],'Color',0.2*[1 1 1]);
line([0.94 0.94],[43.8 41.5],'Color',0.2*[1 1 1]);
line([0.94 0.80],[41.5 41.5],'Color',0.2*[1 1 1]);
line([0.80 0.80],[41.5 43.8],'Color',0.2*[1 1 1]);

line([0.80 0.94],[41.3 41.3],'Color',0.2*[1 1 1]);
line([0.94 0.94],[41.3 39.0],'Color',0.2*[1 1 1]);
line([0.94 0.80],[39.0 39.0],'Color',0.2*[1 1 1]);
line([0.80 0.80],[39.0 41.3],'Color',0.2*[1 1 1]);

line([0.80 0.94],[38.8 38.8],'Color',0.2*[1 1 1]);
line([0.94 0.94],[38.8 36.5],'Color',0.2*[1 1 1]);
line([0.94 0.80],[36.5 36.5],'Color',0.2*[1 1 1]);
line([0.80 0.80],[36.5 38.8],'Color',0.2*[1 1 1]);

% surrounding box
line([0.61 0.95],[36.3 36.3],'Color',0.2*[1 1 1]);
line([0.95 0.95],[36.3 49.0],'Color',0.2*[1 1 1]);
line([0.95 0.61],[49.0 49.0],'Color',0.2*[1 1 1]);
line([0.61 0.61],[49.0 36.3],'Color',0.2*[1 1 1]);

text(0.79,47.6,'Initial WMFT Score','HorizontalAlignment','right');
text(0.79,45.1,'Initial Box-and-Blocks','HorizontalAlignment','right');
text(0.79,42.6,'Max PMTD','HorizontalAlignment','right');
text(0.79,40.1,'Mean PMTD','HorizontalAlignment','right');
text(0.79,37.6,'Initial UEFM Score','HorizontalAlignment','right');

for i = 1:51
    nrow = size(WO_RF,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;
    for j = 1:100
        yy = ymid + ((rand-0.5)*.7); % add noise to y-value
        plot(WO_RF(i,j)/xmax,yy,'.','Color',[205/255 102/255 0],'MarkerSize',5);
    end
%     text(xmin/xmax*8,ymid,rnam4{i},'HorizontalAlignment','left');
    % now plot the box for each value
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',[0 139/255 139/255]);
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255]);
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255]);
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',[0 139/255 139/255]);
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255]); % Median
    
    if i >= 7
        text(Y(i,2)+0.03, ymid, RownamesWO(i),'HorizontalAlignment','left');
    else
        text(Y(i,2)-0.03, ymid, RownamesWO(i),'HorizontalAlignment','right');
    end
end
box off
set(gca,'color','none','YColor','none','Ticklength', [0 0])
xlim([xmin/xmax 1])
xlabel('Normalized Feature Importance')
% title('Variable Importance for Change in Wolf Motor Function')

% plot individual top features
axes('Position',[0.867 0.7 0.054 0.11])
X = [ones(26,1),DF1(:,top5WO(1))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(1)),yWO,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(1)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.87 0.56 0.052 0.11])
X = [ones(26,1),DF1(:,top5WO(2))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(2)),yWO,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(2)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.853 0.42 0.065 0.11])
X = [ones(26,1),DF1(:,top5WO(3))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(3)),yWO,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(3)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.867 0.28 0.0652 0.11])
X = [ones(26,1),DF1(:,top5WO(4))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(4)),yWO,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(4)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.867 0.14 0.0652 0.11])
X = [ones(26,1),DF1(:,top5WO(5))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(5)),yWO,30,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(5)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];
