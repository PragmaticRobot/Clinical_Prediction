clear all
clc
%% This script will do a simple plot of percentage of times a feature was selected by lasso, and also the 

[FM_order_data,Rownames,~] = xlsread('2016dec12_LS_FM_order.csv');
[WO_order_data,~,~] = xlsread('2016dec12_LS_WO_order.csv');

Rownames(1,:) = [];
Rownames(:,2:end) = [];

%% figure out the order for VIs, also normalize VIs

[~,FM_LS_order] = sort(FM_order_data,'descend');
[~,WO_LS_order] = sort(WO_order_data,'descend');

% Switch the features picked 100% of the time so they make more sense
% (initial wolf at top)
WO_LS_order([1 4]) = WO_LS_order([4 1]);
% WO_LS_order([2 4]) = WO_LS_order([4 2]);
WO_LS_order([3 4]) = WO_LS_order([4 3]);

%% housekeeping for top features
[DF1, txt_DF1, ~] = xlsread('2016dec12_df1.csv'); DF1(:,1) = [];
yFM = xlsread('yFM.csv'); yFM(:,1) = [];
yWO = xlsread('yWO.csv'); yWO(:,1) = [];
Names_df1 = txt_DF1(1,:);
top5FM = FM_LS_order(1:5);
top5WO = WO_LS_order(1:5);

%% let's plot
col = 0.7*[1 0 0];
col2 = 0.7*[0 0 1];
figure
clf
axes('Position',[0 0 1 1])
axis off
xx = linspace(15,1,15);
axes('Position',[0.05 0.05 0.45 .95])
hold on;

% create boxes for top features, and a box around all of them
line([0.07 0.47],[4.4 4.4],'Color',0.2*[1 1 1]);
line([0.47 0.47],[4.4 15.2],'Color',0.2*[1 1 1]);
line([0.47 0.07],[15.2 15.2],'Color',0.2*[1 1 1]);
line([0.07 0.07],[15.2 4.4],'Color',0.2*[1 1 1]);

line([0.47 0.93],[15.2 15.0],'Color',0.2*[1 1 1]);
line([0.78 0.47],[10.8 4.4],'Color',0.2*[1 1 1]);

line([0.09 0.21],[13.1 13.1],'Color',0.2*[1 1 1]);
line([0.21 0.21],[13.1 15.0],'Color',0.2*[1 1 1]);
line([0.21 0.09],[15.0 15.0],'Color',0.2*[1 1 1]);
line([0.09 0.09],[15.0 13.1],'Color',0.2*[1 1 1]);

line([0.09 0.21],[10.8 10.8],'Color',0.2*[1 1 1]);
line([0.21 0.21],[10.8 12.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[12.7 12.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[12.7 10.8],'Color',0.2*[1 1 1]);

line([0.09 0.21],[8.8 8.8],'Color',0.2*[1 1 1]);
line([0.21 0.21],[8.8 10.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[10.7 10.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[10.7 8.8],'Color',0.2*[1 1 1]);

line([0.09 0.21],[6.8 6.8],'Color',0.2*[1 1 1]);
line([0.21 0.21],[6.8 8.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[8.7 8.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[8.7 6.8],'Color',0.2*[1 1 1]);

line([0.09 0.21],[4.7 4.7],'Color',0.2*[1 1 1]);
line([0.21 0.21],[4.7 6.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[6.7 6.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[6.7 4.7],'Color',0.2*[1 1 1]);

plot(FM_order_data(FM_LS_order(1:15)),xx,'bd','MarkerSize',7,'MarkerFaceColor','b')
xlim([0 1.1])
ylim([0.5 15.5])
for i = 1:15
    if i == 6 || i == 7 || i == 15 || i == 11
        text(FM_order_data(FM_LS_order(i))+0.02,15-i+1,Rownames(FM_LS_order(i)),'HorizontalAlignment','left');
    else
        text(FM_order_data(FM_LS_order(i))-0.02,15-i+1,Rownames(FM_LS_order(i)),'HorizontalAlignment','right');
    end
end
ax = gca;
yruler = ax.YRuler;
xruler = ax.XRuler;
yruler.Visible = 'off';
xruler.TickLength = 0.2;
set(gca, 'Color', 'None')
% xruler.Axle.VertexData = [0 1.02;0.5 0.5; -1 -1];

text(0.22,14.1,Rownames(FM_LS_order(1)),'HorizontalAlignment','left');
text(0.22,11.5,Rownames(FM_LS_order(2)),'HorizontalAlignment','left');
text(0.22,9.8,Rownames(FM_LS_order(3)),'HorizontalAlignment','left');
text(0.22,7.7,Rownames(FM_LS_order(4)),'HorizontalAlignment','left');
text(0.22,5.5,Rownames(FM_LS_order(5)),'HorizontalAlignment','left');


axes('Position',[0.08 0.85 0.06 0.11])
X = [ones(26,1),DF1(:,top5FM(1))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(1)),yFM,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(1)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.091 0.7 0.04 0.11])
X = [ones(26,1),DF1(:,top5FM(2))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(2)),yFM,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(2)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.09 0.58 0.04 0.11])
X = [ones(26,1),DF1(:,top5FM(3))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(3)),yFM,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(3)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.091 0.45 0.04 0.11])
X = [ones(26,1),DF1(:,top5FM(4))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(4)),yFM,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(4)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.093 0.32 0.04 0.11])
X = [ones(26,1),DF1(:,top5FM(5))]; b = X\yFM; CalcFM = X*b;
scatter(DF1(:,top5FM(5)),yFM,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5FM(5)),CalcFM,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

%% For Wolf now
axes('Position',[0.5 0.05 0.45 .95])
% subplot(1,2,2);
hold on;
plot(WO_order_data(WO_LS_order(1:15)),xx,'bd','MarkerSize',7,'MarkerFaceColor','b')
xlim([0 1.1])
ylim([0.5 15.5])
for i = 1:15
    text(WO_order_data(WO_LS_order(i))-0.02,15-i+1,Rownames(WO_LS_order(i)),'HorizontalAlignment','right');
end
ax = gca;
set(gca, 'Color', 'None')
yruler = ax.YRuler;
xruler = ax.XRuler;
yruler.Visible = 'off';
xruler.TickLength = 0.2;


line([0.07 0.47],[4.4 4.4],'Color',0.2*[1 1 1]);
line([0.47 0.47],[4.4 15.2],'Color',0.2*[1 1 1]);
line([0.47 0.07],[15.2 15.2],'Color',0.2*[1 1 1]);
line([0.07 0.07],[15.2 4.4],'Color',0.2*[1 1 1]);

line([0.47 0.73],[15.2 15.0],'Color',0.2*[1 1 1]);
line([0.83 0.47],[10.8 4.4],'Color',0.2*[1 1 1]);

line([0.09 0.21],[13.1 13.1],'Color',0.2*[1 1 1]);
line([0.21 0.21],[13.1 15.0],'Color',0.2*[1 1 1]);
line([0.21 0.09],[15.0 15.0],'Color',0.2*[1 1 1]);
line([0.09 0.09],[15.0 13.1],'Color',0.2*[1 1 1]);

line([0.09 0.21],[10.8 10.8],'Color',0.2*[1 1 1]);
line([0.21 0.21],[10.8 12.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[12.7 12.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[12.7 10.8],'Color',0.2*[1 1 1]);

line([0.09 0.21],[8.8 8.8],'Color',0.2*[1 1 1]);
line([0.21 0.21],[8.8 10.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[10.7 10.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[10.7 8.8],'Color',0.2*[1 1 1]);

line([0.09 0.21],[6.8 6.8],'Color',0.2*[1 1 1]);
line([0.21 0.21],[6.8 8.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[8.7 8.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[8.7 6.8],'Color',0.2*[1 1 1]);

line([0.09 0.21],[4.7 4.7],'Color',0.2*[1 1 1]);
line([0.21 0.21],[4.7 6.7],'Color',0.2*[1 1 1]);
line([0.21 0.09],[6.7 6.7],'Color',0.2*[1 1 1]);
line([0.09 0.09],[6.7 4.7],'Color',0.2*[1 1 1]);


text(0.22,14.1,Rownames(WO_LS_order(1)),'HorizontalAlignment','left');
text(0.22,11.5,Rownames(WO_LS_order(2)),'HorizontalAlignment','left');
text(0.22,9.8,Rownames(WO_LS_order(3)),'HorizontalAlignment','left');
text(0.22,7.8,Rownames(WO_LS_order(4)),'HorizontalAlignment','left');
text(0.22,5.7,Rownames(WO_LS_order(5)),'HorizontalAlignment','left');


axes('Position',[0.541 0.85 0.04 0.11])
X = [ones(26,1),DF1(:,top5WO(1))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(1)),yWO,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(1)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.542 0.7 0.0652 0.11])
X = [ones(26,1),DF1(:,top5WO(2))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(2)),yWO,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(2)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.54 0.58 0.063 0.11])
X = [ones(26,1),DF1(:,top5WO(3))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(3)),yWO,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(3)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.54 0.45 0.04 0.11])
X = [ones(26,1),DF1(:,top5WO(4))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(4)),yWO,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(4)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];

axes('Position',[0.52 0.32 0.063 0.11])
X = [ones(26,1),DF1(:,top5WO(5))]; b = X\yWO; CalcWO = X*b;
scatter(DF1(:,top5WO(5)),yWO,'MarkerFaceColor',col,'MarkerEdgeColor',col);
hold on
box on
plot(DF1(:,top5WO(5)),CalcWO,'Color',col2,'LineWidth',1)
axis off
ax = gca; ax.XTick = []; ax.YTick = [];


%% beautification

