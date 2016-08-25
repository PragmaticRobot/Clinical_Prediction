% This scripts get all the metrics for all the subjects (not summary, as in
% mean, max and var, but the whole things for every reach
Features = readtable('FeatureSet.csv');
Features(:,1) = [];
subList = {'01','02','03','04','05','06','07','08','09','10','11','12','13',...
    '14','15','16','17','18','19','20','21','22','23','24','25','26','27','28'};
% To get the column names
subMetrics = cell(1,2);
load('D:\Yaz''s Documents\Research\Patient Data\@Mat\AllMetrics_Eval\AllMetrics01.mat')
subMetrics{1,1} = dataset2table(AllMetrics(1:20,:));
subMetrics{1,2} = Features(1,39:end);
FeatNames = cell(34,1);
FeatNames(1:17) = subMetrics{1,1}.Properties.VariableNames;
FeatNames(18:34) = subMetrics{1,2}.Properties.VariableNames;
% FeatNames now contains all 34 column titles for the base feature set
% now we need to put all the data together, we'll just repmat the
% demographic info as many times as the subject center-out reaches
AllMets = zeros(560,34);

for sub = 1:28
    if (sub == 3 || sub == 21)
        continue
    else
        thisSub = subList(sub);
        string = strcat('D:\Yaz''s Documents\Research\Patient Data\@Mat\AllMetrics_Eval\AllMetrics',thisSub,'.mat');
        load(string{1})
        theRows = (1:20)+20*(sub-1);
        AllMets(theRows,1:17) = table2array(dataset2table(AllMetrics(1:20,:)));
    end
end
ToRemove = [41:60,401:420];
AllMets(ToRemove,:)=[];
for sub = 1:26
    theRows = (1:20)+20*(sub-1);
    AllMets(theRows,18:34) = repmat(table2array(Features(sub,39:end)),20,1);
end

save('pqfile.mat','AllMets','FeatNames')

