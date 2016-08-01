% this script plots two trees side by side

%% First Tree Specs
CellText = {'Mean Max Speed','+4.5','Max Initial Movement Ratio',...
    '+0.125','Max Initial Movement Ratio','+4.8','+0.4'};
CellLevel = [1,2,2,3,3,4,4];
CellPos = {'T','T','B','T','B','T','B'};
LvlCutoff = [0.25, 0.86,0.77];
LvlTopSign = {'L','GE','GE'};
LvlCols = {'r','g','b','m'};

%% First Tree Plot
figure(5)
clf
FW = 0.5;
% FW = 3.3;
FH = 1;
% FH = 2.1;
FU = 'normalized';
% FU = 'centimeters';
subplot(1,2,1);
hold on
set(gca,'Position',[0 0 FW FH],'Units',FU)
% set(gca,'Position',[0 0 3.3 2.1],'Units','centimeters')
% first tree limits: 0 ==> 1.65
% second tree limits: 1.65 ==> 3.3

nlvl = max(CellLevel);          % number of levels
boxH = nlvl;                    % height of each box (max)
for i = 1:max(CellLevel)                % Loop over number of levels
    yHmid = FH - (boxH*i)+ 0.5*boxH;    % vetical middle of each box level
    nCell = length(find(CellLevel==i)); % number of cells for that level
    boxW = FW/nCell;                % cell width (max)
    for j = 1:nCell                     % Loop over cells
        yWmid = FW-(boxW*j)+0.5*boxW;   % horizontal center of cell
        plot(yWmid,yHmid,'.','MarkerSize',10,'Color',LvlCols{i})
    end
end
box off
set(gca,'color','none','YColor','none','XColor','none','Ticklength', [0 0])

%% Second Tree Specs
CellText = {'Mean Pre-Movement Speed','Max Number of Speed Peaks','Mass',...
    '+1.44','+5.57','+1','-1'};
CellLevel = [1,2,2,3,3,3,3];
CellPos = {'T','T','B','T','B','T','B'};
LvlCutoff = [0.01, 19.5, 206];
LvlTopSign = {'GE','GE','GE'};
LvlCols = {'r','g','b'};

subplot(1,2,2);
FW = 0.5;
FH = 1;
FU = 'normalized';
set(gca,'Position',[0.5 0 FW FH],'Units',FU)
hold on
nlvl = max(CellLevel);          % number of levels
boxH = FH/max(CellLevel);      % height of each box (max)
for i = 1:nlvl                % Loop over number of levels
    yHmid = FH - (boxH*i)+ 0.5*boxH;    % vetical middle of each box level
    nCell = length(find(CellLevel==i)); % number of cells for that level
    boxW = FW/nCell;                % cell width (max)
    for j = 1:nCell                     % Loop over cells
        yWmid = FW-(boxW*j)+0.5*boxW;   % horizontal center of cell
        plot(yWmid,yHmid,'.','MarkerSize',10,'Color',LvlCols{i})
    end
end
box off
set(gca,'color','none','YColor','none','XColor','none','Ticklength', [0 0])
% 
% ymid = (nrow-i)+0.5;
% xx = [xmin xmax xmax xmin];
% yy = [ymid-0.5 ymid-0.5 ymid+0.5 ymid+0.5];
% col = map(countFM(i),:);
% patch(xx,yy,col,'EdgeColor','none');

%% Try again
% figure(7)
% clf
% set(gca,'Position',[0 0 .5 1])
% 
% % level 1
% xx=[35 65 65 35];
% yy=[80 80 90 90];
% col = [0 139/255 139/255];
% patch(xx,yy,col,'EdgeColor','none')
% 
% % level 2 - 2 cells
% xx=[65 85 85 65];
% yy=[60 60 70 70];
% patch(xx,yy,col,'EdgeColor','none')
% 
% xx=[30 40 40 30];
% yy=[60 60 70 70];
% patch(xx,yy,col,'EdgeColor','none')
% 
% xlim([0 100])
% ylim([0 100])
% level 3 - 2 cells
% xx=[0.05 .15 0.15 0.05];
% yy=[0.35 0.35 0.45 0.45];
% patch(xx,yy,col,'EdgeColor','none')
% 
% xx=[0.175 .2 0.275 0.1];
% yy=[0.35 0.35 0.45 0.45];
% patch(xx,yy,col,'EdgeColor','none')
% 
% xx=[0.1 .2 0.2 0.1];
% yy=[0.5 0.5 0.6 0.6];
% patch(xx,yy,col,'EdgeColor','none')
% 
% xx=[0.1 .2 0.2 0.1];
% yy=[0.5 0.5 0.6 0.6];
% patch(xx,yy,col,'EdgeColor','none')