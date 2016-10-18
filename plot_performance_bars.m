function [ ] = plot_performance_bars(Means, Stds, locs, ylcode,plotCode,outCode)
%plot_performance_bars plots bar graphs to compare performance metrics of
%different models. The metrics are represented in ylcode. The type of
%comparison is represented by plotCode. locs is a vector with the locations
%of the bars
% ylcode legend:
% 1: RMSE / 2: R2 / 3: Adjusted R2 / 4: Slope
% plotCode legend:
% 1: Comparing the different approaches to lasso (v8 - v13)
% 2: Comparing lasso and RF, different fold approaches (v11-v13 lasso, v8-v10 RF)
% 3: Comparing picked lasso and RF (vXX lasso and vXX RF), Fugl-Meyer and Wolf

switch plotCode
    case 1
        % when comparing the different lasso models, we plot 12 bars, 6 for
        % each type of cross validation (best and worst performance.
        colors = repmat([0.933 0.172 0.172;0 0.545 0.545],[length(locs)/2,1]);

        figure
        hold on

        for i = 1:length(locs)
            bar(locs(i),Means(i),0.5,'FaceColor',colors(i,:));
            plot( [locs(i) locs(i)], Means(i)+Stds(i)*[1 -1], 'color','k',  'lineWidth',1.5);
        end
        
        ax = gca;
        ax.XTick = [4.5 15.5];
        ax.XTickLabel = {'Externally controlled CV','Best Performance CV'};
        set(ax,'FontSize',18)
        textlocs = [1,4,7,12,15,18];
        versions = {'4-fold','Leave-One-Out','13-fold','4-fold','Leave-One-Out','13-fold'};
        
        switch ylcode
            case 1
                ylabel('RMSE')
                if outCode == 2 % if it's RMSE and Wolf
                    for i = 1:length(textlocs)
                        text(textlocs(i),6.6,versions(i),'FontSize',15);
                    end
                    ylim([0 7])
                    set(ax,'YTick',0:0.5:7);
                else
                    for i = 1:length(textlocs)
                        text(textlocs(i),5.5,versions(i),'FontSize',15);
                    end
                    ylim([0 6.5])
                    set(ax,'YTick',0:0.5:5);
                end
            case 2
                ylabel('R^2')
                for i = 1:length(textlocs)
                    text(textlocs(i),1,versions(i),'FontSize',15);
                end
                ylim([0 1.2])
                set(ax,'YTick',0:0.1:1);
            case 3
                ylabel('Adjusted R^2')
                for i = 1:length(textlocs)
                    text(textlocs(i),1,versions(i),'FontSize',15);
                end
                ylim([0 1.2])
                set(ax,'YTick',0:0.1:1);
            case 4
                ylabel('Slope')
                for i = 1:length(textlocs)
                    text(textlocs(i),1,versions(i),'FontSize',15);
                end
                ylim([0 1.2])
                set(ax,'YTick',0:0.1:1);
        end
        
        textlabel = repmat({'First order lasso';'Second order lasso'},[6,1]);
        
        for i = 1:length(locs)
            text(locs(i)-0.02,0.03,textlabel(i),'Rotation',90,'FontSize',15)
        end
        
    case 2
        colors = repmat([0.933 0.172 0.172;0 0.545 0.545;1 0.5 0;0.4314 0.5961 0.5294],[length(locs)/4,1]);

        figure
        hold on

        for i = 1:length(locs)
            bar(locs(i),Means(i),0.5,'FaceColor',colors(i,:));
            plot( [locs(i) locs(i)], Means(i)+Stds(i)*[1 -1], 'color','k',  'lineWidth',1.5);
        end
        
        ax = gca;
        ax.XTick = [];
        ax.XTickLabel = [];
        set(ax,'FontSize',18)
        textlocs = [2,7,12];
        versions = {'4-fold','Leave-One-Out','13-fold'};
        
        switch ylcode
            case 1
                ylabel('RMSE')
                if outCode == 2
                    for i = 1:length(textlocs)
                        text(textlocs(i),6.6,versions(i),'FontSize',15);
                    end
                    ylim([0 7])
                    set(ax,'YTick',0:0.5:7);
                else
                    for i = 1:length(textlocs)
                        text(textlocs(i),5.5,versions(i),'FontSize',15);
                    end
                    ylim([0 6.5])
                    set(ax,'YTick',0:0.5:5);
                end
            case 2
                ylabel('R^2')
                for i = 1:length(textlocs)
                    text(textlocs(i),1,versions(i),'FontSize',15);
                end
                ylim([0 1.2])
                set(ax,'YTick',0:0.1:1);
            case 3
                ylabel('Adjusted R^2')
                for i = 1:length(textlocs)
                    text(textlocs(i),1,versions(i),'FontSize',15);
                end
                ylim([0 1.2])
                set(ax,'YTick',0:0.1:1);
            case 4
                ylabel('Slope')
                for i = 1:length(textlocs)
                    text(textlocs(i),1,versions(i),'FontSize',15);
                end
                ylim([0 1.2])
                set(ax,'YTick',0:0.1:1);
        end
        
        textlabel = repmat({'First order lasso';'Second order lasso';'First order random forests';'Second order random forests'},[3,1]);
        
        for i = 1:length(locs)
            text(locs(i)-0.02,0.03,textlabel(i),'Rotation',90,'FontSize',15)
        end
        
    otherwise
end
% legend([p1 p2 p3 p4],{'First order lasso';'Second order lasso';...
%     'First order random forests';'Second order random forests'},'FontSize',12,'Location','north','Box','off')


end

