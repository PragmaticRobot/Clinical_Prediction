clear all
clc

%% load the VI matrices

FM_VI_v8 = csvread('FM_VI_v8.csv',1,1);
FM_VI_v9 = csvread('FM_VI_v9.csv',1,1);
FM_VI_v10 = csvread('FM_VI_v10.csv',1,1);

[~,Titles,~] = xlsread('FM_VI_v8.csv');
Titles(1,:) = [];
Titles(:,2:end) = [];

%% 8 features per sheet

Plot_VIs(FM_VI_v8(1:8,:),Titles(1:8))
Plot_VIs(FM_VI_v8(9:16,:),Titles(9:16))
Plot_VIs(FM_VI_v8(17:24,:),Titles(17:24))
Plot_VIs(FM_VI_v8(25:32,:),Titles(25:32))
Plot_VIs(FM_VI_v8(33:40,:),Titles(33:40))
Plot_VIs(FM_VI_v8(41:48,:),Titles(41:48))
Plot_VIs(FM_VI_v8(49:54,:),Titles(49:54))

Plot_VIs(FM_VI_v9(1:8,:),Titles(1:8))
Plot_VIs(FM_VI_v9(9:16,:),Titles(9:16))
Plot_VIs(FM_VI_v9(17:24,:),Titles(17:24))
Plot_VIs(FM_VI_v9(25:32,:),Titles(25:32))
Plot_VIs(FM_VI_v9(33:40,:),Titles(33:40))
Plot_VIs(FM_VI_v9(41:48,:),Titles(41:48))
Plot_VIs(FM_VI_v9(49:54,:),Titles(49:54))

Plot_VIs(FM_VI_v10(1:8,:),Titles(1:8))
Plot_VIs(FM_VI_v10(9:16,:),Titles(9:16))
Plot_VIs(FM_VI_v10(17:24,:),Titles(17:24))
Plot_VIs(FM_VI_v10(25:32,:),Titles(25:32))
Plot_VIs(FM_VI_v10(33:40,:),Titles(33:40))
Plot_VIs(FM_VI_v10(41:48,:),Titles(41:48))
Plot_VIs(FM_VI_v10(49:54,:),Titles(49:54))
