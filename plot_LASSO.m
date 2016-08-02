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
%% plotting LASSO Fugl-Meyer

% box limits
Y = prctile(FM_LASSO,[25 50 75],2); % calculate 25% 50% 75% for each feature

xmax = ceil(max(max(FM_LASSO)));
xmin = floor(min(min(FM_LASSO)));
figure(1)
clf
set(gca,'Position',[.05 .1 .85 .75])
hold on
for i = 1:size(FM_LASSO,1)
    nrow = size(FM_LASSO,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;   
    xx = [xmin xmax xmax xmin];
    yy = [ymid-0.5 ymid-0.5 ymid+0.5 ymid+0.5];
    col = map(countFM(i),:);
    patch(xx,yy,col,'EdgeColor','none');
    for j = 1:100
        if isnan(FM_LASSO(i,j))
            continue
        end
        yy = ymid + (rand-0.5); % add noise to y-value
        plot(FM_LASSO(i,j),yy,'.','Color',mapComp(countFM(i),:),'MarkerSize',5);
    end
    text(xmax-0.3,ymid,rnam1{i},'HorizontalAlignment','right');
    % now plot the box for each value
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',mapComp(countFM(i),:));
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',mapComp(countFM(i),:));
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',mapComp(countFM(i),:));
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',mapComp(countFM(i),:));
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',mapComp(countFM(i),:)); % Median
    h1 = plot([0 0],[0 nrow],'Color',[0.016 0.184 0.184],'LineStyle','--');
    h1.Color(4) = 1;
end
colormap(map)
c1 = colorbar('Ticks',[0,.20,.40,.60,.80,1.00],'TickLabels',{'0','20','40',...
    '60','80','100'},'Location','east','Position',[0.92 0.13 0.02 0.7]);
c1.Label.String = 'Selection Frequency';
box off
set(gca,'YColor','none','Color','none')
xlim([xmin xmax])
title('LASSO Features for Predicting Change in UEFM')
xlabel('Coefficient')
axis tight

%% Plotting LASSO Wolf
% box limits
Y = prctile(WO_LASSO,[25 50 75],2); % calculate 25% 50% 75% for each feature

xmax = ceil(max(max(WO_LASSO)));
xmin = floor(min(min(WO_LASSO)));
figure(2)
clf
set(gca,'Position',[.05 .1 .85 .75])
hold on
for i = 1:size(WO_LASSO,1)
    nrow = size(WO_LASSO,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;   
    xx = [xmin xmax xmax xmin];
    yy = [ymid-0.5 ymid-0.5 ymid+0.5 ymid+0.5];
    col = map2(countWO(i),:);
    patch(xx,yy,col,'EdgeColor','none');
    for j = 1:100
        if isnan(WO_LASSO(i,j))
            continue
        end
        yy = ymid + (rand-0.5); % add noise to y-value
        plot(WO_LASSO(i,j),yy,'.','Color',map2Comp(countWO(i),:),'MarkerSize',5);
    end
    text(xmax-1,ymid,rnam3{i},'HorizontalAlignment','right');
    % now plot the box for each value
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',map2Comp(countWO(i),:));
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',map2Comp(countWO(i),:));
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',map2Comp(countWO(i),:));
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',map2Comp(countWO(i),:));
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',map2Comp(countWO(i),:)); % Median
    h1 = plot([0 0],[0 nrow],'Color',[0.016 0.184 0.184],'LineStyle','--');
    h1.Color(4) = 1;
end
colormap(map2)
c1 = colorbar('Ticks',[0,.20,.40,.60,.80,1.00],'TickLabels',{'0','20','40',...
    '60','80','100'},'Location','east','Position',[0.92 0.13 0.02 0.7]);
c1.Label.String = 'Selection Frequency';
box off
set(gca,'YColor','none','Color','none')
xlim([xmin xmax])
title('LASSO Features for Predicting Change in Wolf Motor Function')
xlabel('Coefficient')
axis tight