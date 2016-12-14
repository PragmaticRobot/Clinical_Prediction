
%% This script reads the output from R and then builds bottom-up models to compare the different kinds of models

clear all
clc
%% Load the basic data
[DF1, txt_DF1, ~] = xlsread('df1.csv');     % load the base feature set
% Remember features 43:51 are categorical
yFM = xlsread('yFM.csv');                   % load delta FM
yWO = xlsread('yWO.csv');                   % load delta wolf
Names_df1 = txt_DF1(1,:);                   % feature names
DF1(:,1) = [];                              % remove extra columns
yFM(:,1) = [];                              % same
yWO(:,1) = [];                              % same

% the top lasso features from Fugl-Meyer, to be used for the curve fitting
% toolbox
MeanMaxSpd      = DF1(:,5);
VarSpdRto       = DF1(:,33);
MonthsPost      = DF1(:,42);
Age             = DF1(:,39);
MaxMaxSpd       = DF1(:,18);

%% Find the feature order and set up parameters

nSubj   = 26;                 % how many subjects in your dataset?
nFeats  = 20;                 % max number of features on x-axis of final plot
x       = 1;                  % how many do you want to leave out?
xx      = 1:1:26;             % this will be trashed soon
Folds   = nchoosek(xx,x);     % all the possible combination of Leave x out
nReps   = size(Folds,1);      % make sure to cover all possibilities
nFold   = 26;                 % nFold cross validation for lasso
clear xx                      % remove the trash

LS_FM_order = xlsread('LS_FM_order.csv');       % winning features for FM using lasso
LS_WO_order = xlsread('LS_WO_order.csv');       % winning features for WO using lasso
[~,FM_LS_ind] = sort(LS_FM_order,'descend');
[~,WO_LS_ind] = sort(LS_WO_order,'descend');

%% Set up placeholders for things
ThisMat_FM_LS = zeros(nSubj,nFeats);

FM_LS_L = zeros(nReps,nFeats);
FM_LS_Q = zeros(nReps,nFeats);

ThisBaseQuad = zeros(nSubj,nFeats);
IsCateg = zeros(nFeats,1);

for i = 1:nFeats
    if FM_LS_ind(i) >=43 && FM_LS_ind(i) <= 51
        IsCateg(i) = 1;
    else
        IsCateg(i) = 0;
    end
end

%% Build the models and use

thisOut = yFM;              % what's our output metric here?
for i = 1:nFeats
    ThisMat_FM_LS(:,i) = DF1(:,FM_LS_ind(i));
    for j = 1:nReps
        ThisCV = ThisMat_FM_LS(:,1:i);
        %         ThisCV(Folds(j,:),:) = [];
        %         thisOut(Folds(j,:),:) = []; % features 1:i, remove jth subject
        [trash,trash2] = lasso(ThisCV, thisOut, 'Alpha',1,'CV',nFold);	  % fit lasso model with nFold cross-validation
        % trash contains the p by L regression coefficients, trash2
        % contains fit information, which will be used now to get the
        % winning model
        Win_mod = trash(:,trash2.Index1SE);
        xx = ThisCV*Win_mod; % outcome values for winning model
        trash3 = fitlm(xx,thisOut,'linear');
        FM_LS_L(j,i) = trash3.Rsquared.Ordinary;
        %% For quad feature set
        % set up the jth column first, if it's categorical
        if FM_LS_ind(i) >= 43 && FM_LS_ind(i) <= 51
            ThisBaseQuad(:,i) = ones(nSubj,1)-ThisCV(:,i);
        else
            ThisBaseQuad(:,i) = ThisCV(:,i);
        end
        CleanBase = ThisBaseQuad(:,1:i);
        if i >= 3 && i < 5
            ThisQuadCV = x2fx(CleanBase,'quadratic',3);
        elseif i >= 5 && i < 15
            ThisQuadCV = x2fx(CleanBase,'quadratic',[3,5]);
        elseif i >= 15 && i <= 20
            ThisQuadCV = x2fx(CleanBase,'quadratic',[3,5,15]);
            % Make sure you add any extra categorical cases if nFeat > 20
        else
            ThisQuadCV = x2fx(CleanBase,'quadratic');
        end
        [trash,trash2] = lasso(ThisQuadCV, thisOut, 'Alpha',1,'CV',nFold);	  % fit lasso model with nFold cross-validation
        % trash contains the p by L regression coefficients, trash2
        % contains fit information, which will be used now to get the
        % winning model
        Win_mod = trash(:,trash2.Index1SE);
        xx = ThisQuadCV*Win_mod; % outcome values for winning model
        trash3 = fitlm(xx,thisOut,'linear');
        FM_LS_Q(j,i) = trash3.Rsquared.Ordinary;
    end
end

%% Now let's plot this
colors(1,:) = [0.933 0.172 0.172]; colors(2,:) = [0 0.545 0.545]; colors(3,:) = [1 0.5 0]; colors(4,:) = [110, 152, 135]/255;

xx2 = 1:20;
figure
clf
hold on

xx = repmat(1:20,26,1);
for i = 1:20;
    scatter(xx(:,i),FM_LS_L(:,i),1,colors(1,:))
    scatter(xx(:,i),FM_LS_Q(:,i),1,colors(2,:))
    scatter(xx(:,i),FM_RF_L(:,i),1,colors(3,:))
    scatter(xx(:,i),FM_RF_Q(:,i),1,colors(4,:))
end

q1 = plot(xx2,mean(FM_LS_L,1),'Color',colors(1,:),'LineWidth',1.5);
q2 = plot(xx2,mean(FM_LS_Q,1),'Color',colors(2,:),'LineWidth',1.5);
q3 = plot(xx2,mean(FM_RF_L,1),'Color',colors(3,:),'LineWidth',1.5);
q4 = plot(xx2,mean(FM_RF_Q,1),'Color',colors(4,:),'LineWidth',1.5);

legend([q1;q2;q3;q4],{'1st order lasso';'2nd order lasso';'1st order RF';'2nd order RF'})