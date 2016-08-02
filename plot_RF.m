% This script takes the results of the R models (LASSO, RF) with ranks for
% both Fugl-Meyer (FM) and Wolf (WO).

%% data prep
load('R_results.mat')

% convert data type so it's more usable

gg = struct2cell(FM_Lin_LASSO_sort);
gg = gg';
FM_LASSO = cell2mat(gg);

gg = struct2cell(FM_Lin_RF_sort);
gg = gg';
FM_RF = cell2mat(gg);

gg = struct2cell(Wolf_Lin_LASSO_sort);
gg = gg';
WO_LASSO = cell2mat(gg);
trash = WO_LASSO(2,:);
WO_LASSO(2,:) = WO_LASSO(11,:);
WO_LASSO(11,:) = trash;

trash = rnam3{2};
rnam3{2} = rnam3{11};
rnam3{11} = trash;

gg = struct2cell(Wolf_Lin_RF_sort);
gg = gg';
WO_RF = cell2mat(gg);

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
%% plotting Random Forests Fugl-Meyer

% box limits
Y = prctile(FM_RF,[25 50 75],2); % calculate 25% 50% 75% for each feature
xmax = ceil(max(max(FM_RF)));
xmin = floor(min(min(FM_RF(1:15,:))));
Y = Y/xmax;
figure(1)
clf
subplot(1,2,1);
set(gca,'Position',[0.05 .1 .45 .8])
hold on
for i = 1:15
    nrow = size(FM_RF,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;
    for j = 1:100
        yy = ymid + ((rand-0.5)*.7); % add noise to y-value
        plot(FM_RF(i,j)/xmax,yy,'.','Color',[0 139/255 139/255],'MarkerSize',5);
    end
    text(xmin/xmax*2.1,ymid,rnam2{i},'HorizontalAlignment','left');
    % now plot the box for each value
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',[0 139/255 139/255]);
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255]);
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255]);
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',[0 139/255 139/255]);
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255]); % Median
end
box off
set(gca,'color','none','YColor','none','Ticklength', [0 0])
xlim([xmin/xmax 1])
xlabel('Normalized Feature Importance')
title('Variable Importance for Change in Fugl-Meyer')
%% plotting Random Forests Wolf

% box limits
Y = prctile(WO_RF,[25 50 75],2); % calculate 25% 50% 75% for each feature
xmax = ceil(max(max(WO_RF)));
xmin = floor(min(min(WO_RF(1:15,:))));
Y = Y/xmax;
% figure(2)
% clf
subplot(1,2,2);
set(gca,'Position',[.55 .1 .4 .8])
hold on
for i = 1:15
    nrow = size(WO_RF,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;
    for j = 1:100
        yy = ymid + ((rand-0.5)*.7); % add noise to y-value
        plot(WO_RF(i,j)/xmax,yy,'.','Color',[205/255 102/255 0],'MarkerSize',5);
    end
    text(xmin/xmax*8,ymid,rnam2{i},'HorizontalAlignment','left');
    % now plot the box for each value
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',[205/255 102/255 0]);
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0]);
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0]);
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',[205/255 102/255 0]);
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0]); % Median
end
box off
set(gca,'color','none','YColor','none','Ticklength', [0 0])
xlim([xmin/xmax 1])
xlabel('Normalized Feature Importance')
title('Variable Importance for Change in Wolf Motor Function')
