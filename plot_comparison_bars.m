function [ ] = plot_comparison_bars( Means, Stds,yl)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

colors = zeros(4,3);
colors(1,:) = [1 0 0]; colors(2,:) = [0 0 1]; colors(3,:) = [1 0.7 0]; colors(4,:) = [0 1 0];

% colors(1,:) = [0.933 0.172 0.172]; colors(2,:) = [0 0.545 0.545]; colors(3,:) = [1 0.5 0]; colors(4,:) = [110, 152, 135]/255;

figure
hold on
p1 = bar(1,Means(1),0.5,'FaceColor',colors(1,:));
p2 = bar(2,Means(2),0.5,'FaceColor',colors(2,:));
p3 = bar(3,Means(3),0.5,'FaceColor',colors(3,:));
p4 = bar(4,Means(4),0.5,'FaceColor',colors(4,:));

plot( [1 1], Means(1)+Stds(1)*[1 -1], 'color','k',  'lineWidth',1);
plot( [2 2], Means(2)+Stds(2)*[1 -1], 'color','k',  'lineWidth',1);
plot( [3 3], Means(3)+Stds(3)*[1 -1], 'color','k',  'lineWidth',1);
plot( [4 4], Means(4)+Stds(4)*[1 -1], 'color','k',  'lineWidth',1);

bar(6,Means(5),0.5,'FaceColor',colors(1,:));
bar(7,Means(6),0.5,'FaceColor',colors(2,:));
bar(8,Means(7),0.5,'FaceColor',colors(3,:));
bar(9,Means(8),0.5,'FaceColor',colors(4,:));

plot( [6 6], Means(5)+Stds(5)*[1 -1], 'color','k',  'lineWidth',1);
plot( [7 7], Means(6)+Stds(6)*[1 -1], 'color','k',  'lineWidth',1);
plot( [8 8], Means(7)+Stds(7)*[1 -1], 'color','k',  'lineWidth',1);
plot( [9 9], Means(8)+Stds(8)*[1 -1], 'color','k',  'lineWidth',1);

ax = gca;
ax.XTick = [2.5 7.5];
ax.XTickLabel = {'Fugl-Meyer','Wolf Motor Function'};
set(ax,'FontSize',18)
ylabel(yl)
text(0.98,0.02,'First order lasso','Rotation',90,'FontSize',15)
text(1.98,0.02,'Second order lasso','Rotation',90,'FontSize',15,'Color','w')
text(2.98,0.02,'First order random forests','Rotation',90,'FontSize',15)
text(3.98,0.02,'Second order random forests','Rotation',90,'FontSize',15)

text(5.98,0.02,'First order lasso','Rotation',90,'FontSize',15)
text(6.98,0.02,'Second order lasso','Rotation',90,'FontSize',15,'Color','w')
text(7.98,0.02,'First order random forests','Rotation',90,'FontSize',15)
text(8.98,0.02,'Second order forests','Rotation',90,'FontSize',15)

% legend([p1 p2 p3 p4],{'First order lasso';'Second order lasso';...
%     'First order random forests';'Second order random forests'},'FontSize',12,'Location','north','Box','off')

end

