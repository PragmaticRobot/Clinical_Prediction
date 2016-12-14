clear all
clc

Features = readtable('FeatureSet.csv');
Features(:,1) = [];

affected = str2num(cell2mat(Features.affected));
same = str2num(cell2mat(Features.same));

XF = Features(:,[39,40,33,5,42,9,18]);

corrplot(XF,'type','Spearman','varNames',XF.Properties.VariableNames,'testR','on')

XW = Features(:,[54,19,35,5,55,24,12,53]);

corrplot(XW,'type','Spearman','varNames',XW.Properties.VariableNames,'testR','on')