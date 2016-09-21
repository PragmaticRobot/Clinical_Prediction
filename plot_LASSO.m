% This script takes the results of the R models (LASSO, RF) with ranks for
% both Fugl-Meyer (FM) and Wolf (WO).
clear all
clc
%% data prep
load('R_results_4FoldSep16.mat')
load('countsF_4FoldSep16.mat')
Features = readtable('FeatureSet.csv');
Features(:,1) = [];

FM_count = countsFM.Count;
WO_count = countsWO.Count;

[BFM,IFM] = sort(FM_count,'descend');
[BWO,IWO] = sort(WO_count,'descend');

IFM = IFM(1:21);
IWO = IWO(1:1);

FM_Feats = Features(:,IFM);
WO_Feats = Features(:,IWO);

meansFM = mean(table2array(FM_Feats));
meansWO = mean(table2array(WO_Feats));
% convert data type so it's more usable

gg = struct2cell(FM_Lin_LASSO_sort);
gg = gg';
FM_LASSO = cell2mat(gg);
trash = repmat(meansFM',1,100);
FM_LASSO = FM_LASSO.*trash;

gg = struct2cell(FM_Lin_RF_sort);
gg = gg';
FM_RF = cell2mat(gg);

gg = struct2cell(Wolf_Lin_LASSO_sort);
gg = gg';
WO_LASSO = cell2mat(gg);
% trash = WO_LASSO(2,:);
% WO_LASSO(2,:) = WO_LASSO(11,:);
% WO_LASSO(11,:) = trash;
% trash = meansWO(2);
% meansWO(2) = meansWO(11);
% meansWO(11) = trash;
trash = repmat(meansWO',1,100);
WO_LASSO = WO_LASSO.*trash;

% trash = rnam3{2};
% rnam3{2} = rnam3{11};
% rnam3{11} = trash;

gg = struct2cell(Wolf_Lin_RF_sort);
gg = gg';
WO_RF = cell2mat(gg);
FMLS = FM_LASSO;
WOLS = WO_LASSO;
% do logs for FM_LASSO and WO_LASSO
for i = 1:size(FM_LASSO,1)
    for j = 1:size(FM_LASSO,2)
        if FM_LASSO(i,j) > 1
            FM_LASSO(i,j) = log(FM_LASSO(i,j));
        elseif FM_LASSO(i,j) < 1 && FM_LASSO(i,j) > 0
            FM_LASSO(i,j) = -log(FM_LASSO(i,j));
        elseif FM_LASSO(i,j) < 0 && FM_LASSO(i,j) > -1
            FM_LASSO(i,j) = log(abs(FM_LASSO(i,j)));
        elseif FM_LASSO(i,j) < -1
            FM_LASSO(i,j) = -log(abs(FM_LASSO(i,j)));
        elseif isnan(FM_LASSO(i,j))
            FM_LASSO(i,j) = NaN;
        end
    end
end

for i = 1:size(WO_LASSO,1)
    for j = 1:size(WO_LASSO,2)
        if WO_LASSO(i,j) > 1
            WO_LASSO(i,j) = log(WO_LASSO(i,j));
        elseif WO_LASSO(i,j) < 1 && WO_LASSO(i,j) > 0
            WO_LASSO(i,j) = -log(WO_LASSO(i,j));
        elseif WO_LASSO(i,j) < 0 && WO_LASSO(i,j) > -1
            WO_LASSO(i,j) = log(abs(WO_LASSO(i,j)));
        elseif WO_LASSO(i,j) < -1
            WO_LASSO(i,j) = -log(abs(WO_LASSO(i,j)));
        elseif isnan(WO_LASSO(i,j))
            WO_LASSO(i,j) = NaN;
        end
    end
end
% Use name abbreviations:
% rnam1{6} = 'Max NSP (-)';
% rnam1{7} = 'Init BB (+)';
% rnam1{11} = 'Mean PMS';
% rnam1{12} = 'Mean IDE';
% rnam1{13} = 'Max MAPR';
% rnam1{14} = 'Var MAPR';
% rnam1{15} = 'EA';
% rnam1{16} = 'Mean PMTD';
% rnam1{8} = 'Months Post-Stroke (-)';
% 
% rnam3{1} = 'Init WO';
% rnam3{6} = 'Age (+)';
% rnam3{7} = 'Var NSP (+)';
% rnam3{8} = 'Max PLR (+)';
% rnam3{9} = 'Max HPL';
% rnam3{10} = 'Mean MAPR';
% rnam3{12} = 'Var MAPR';
% rnam3{16} = 'Max IDE';
% rnam3{17} = 'Hemorrhagic Stroke (+)';
% rnam3{19} = 'Mean PMTD';
% rnam3{20} = 'EA';
% rnam3{22} = 'Var PMS';
% rnam3{23} = 'Init BB';

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
ax1 = subplot(1,2,1);
colormap(ax1,map)
c1 = colorbar('Ticks',[0,.20,.40,.60,.80,1.00],'TickLabels',{'0','20','40',...
    '60','80','100'},'Location','east','Position',[0.46 0.13 0.01 0.7]);
c1.Label.String = 'Selection Frequency';
box off
set(gca,'YColor','none','Color','none')
xlim([xmin xmax])
title('LASSO Features for Predicting Change in UEFM')
xlabel('Coefficient')
axis tight
set(gca,'Position',[.05 .1 .4 .75])
LL = xmax - xmin;
hold on
for i = 1:size(FM_LASSO,1)
    nrow = size(FM_LASSO,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;   
    xmin = xmax-((countFM(i)*LL)/100);
    xx = [xmin xmax xmax xmin];
    yy = [ymid-0.5 ymid-0.5 ymid+0.5 ymid+0.5];
    col = map(countFM(i),:);
    patch(xx,yy,col,'EdgeColor','none');
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',[205/255 102/255 0],'LineWidth',2);
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0],'LineWidth',2);
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0],'LineWidth',2);
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',[205/255 102/255 0],'LineWidth',2);
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',[205/255 102/255 0],'LineWidth',2); % Median
    for j = 1:100
        if isnan(FM_LASSO(i,j))
            continue
        end
        yy = ymid + (rand-0.5); % add noise to y-value
        plot(FM_LASSO(i,j),yy,'.','Color',mapComp(countFM(i),:),'MarkerSize',5);
    end
    text(xmax-0.3,ymid,rnam1{i},'HorizontalAlignment','right','FontWeight','bold');
    % now plot the box for each value
    h1 = plot([0 0],[0 nrow],'Color',[0.016 0.184 0.184],'LineStyle','--');
    h1.Color(4) = 1;
end



%% Plotting LASSO Wolf
% box limits
Y = prctile(WO_LASSO,[25 50 75],2); % calculate 25% 50% 75% for each feature

xmax = ceil(max(max(WO_LASSO)));
xmin = floor(min(min(WO_LASSO)));
% figure(1)
% clf
ax2 = subplot(1,2,2);
colormap(ax2,map2)
c2 = colorbar('Ticks',[0,.20,.40,.60,.80,1.00],'TickLabels',{'0','20','40',...
    '60','80','100'},'Location','east','Position',[0.96 0.13 0.01 0.7]);
c2.Label.String = 'Selection Frequency';
box off
set(gca,'YColor','none','Color','none')
xlim([xmin xmax])
title('LASSO Features for Predicting Change in Wolf Motor Function')
xlabel('Coefficient')
axis tight
set(gca,'Position',[.52 .1 .43 .75])
LL = xmax - xmin;
hold on
for i = 1:size(WO_LASSO,1)
    nrow = size(WO_LASSO,1);
    % center of box at nrow-i+0.5
    ymid = (nrow-i)+0.5;  
    xmin = xmax-((countWO(i)*LL)/100);
    xx = [xmin xmax xmax xmin];
    yy = [ymid-0.5 ymid-0.5 ymid+0.5 ymid+0.5];
    col = map2(countWO(i),:);
    patch(xx,yy,col,'EdgeColor','none');
    line([Y(i,1) Y(i,3)],[ymid-0.25 ymid-0.25],'Color',[0 139/255 139/255],'LineWidth',2);
    line([Y(i,1) Y(i,1)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255],'LineWidth',2);
    line([Y(i,3) Y(i,3)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255],'LineWidth',2);
    line([Y(i,3) Y(i,1)],[ymid+0.25 ymid+0.25],'Color',[0 139/255 139/255],'LineWidth',2);
    line([Y(i,2) Y(i,2)],[ymid-0.25 ymid+0.25],'Color',[0 139/255 139/255],'LineWidth',2); % Median
    for j = 1:100
        if isnan(WO_LASSO(i,j))
            continue
        end
        yy = ymid + (rand-0.5); % add noise to y-value
        plot(WO_LASSO(i,j),yy,'.','Color',map2Comp(countWO(i),:),'MarkerSize',5);
    end
    text(xmax-2,ymid,rnam3{i},'HorizontalAlignment','right','FontWeight','bold');
    % now plot the box for each value
    h1 = plot([0 0],[0 nrow],'Color',[0.016 0.184 0.184],'LineStyle','--');
    h1.Color(4) = 1;
end
