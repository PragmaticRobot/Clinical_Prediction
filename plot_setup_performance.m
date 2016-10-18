clear all
clc

load('performance_metrics.mat')

% Info to send to the plotting function:
% bar locations
% number of different types of bars (lasso & RF or just lasso)
% Means
% standard deviations (wings)
% bar locations
% code for which output it is (RMSE, R2, AR2, Slope)
% code for the type of comparison
% code for which outcome

% first let's plot all the lasso types (v8 - v13)
% locs will leave a gap of 1 between each version outcome
% with a bigger gap between outside-controlled cv and proper cv

locs = [1,2,4,5,7,8,12,13,15,16,18,19];
means = RMSE_M(1:12,1);
stds = RMSE_STD(1:12,1);
plot_performance_bars(means, stds,locs,1,1,1)

locs = [1,2,4,5,7,8,12,13,15,16,18,19];
means = R2_M(1:12,1);
stds = R2_STD(1:12,1);
plot_performance_bars(means, stds,locs,2,1,1)

% Now let's compare best lasso and RF performance in both RMSE and R2
locs = [1,2,3,4,6,7,8,9,11,12,13,14];
means = RMSE_M([7,8,13,14,9,10,15,16,11,12,17,18],1);
stds = RMSE_STD([7,8,13,14,9,10,15,16,11,12,17,18],1);
plot_performance_bars(means, stds,locs,1,2,1)

locs = [1,2,3,4,6,7,8,9,11,12,13,14];
means = R2_M([7,8,13,14,9,10,15,16,11,12,17,18],1);
stds = R2_STD([7,8,13,14,9,10,15,16,11,12,17,18],1);
plot_performance_bars(means, stds,locs,2,2,1)

%% Do the same for Wolf
locs = [1,2,4,5,7,8,12,13,15,16,18,19];
means = RMSE_M(1:12,2);
stds = RMSE_STD(1:12,2);
plot_performance_bars(means, stds,locs,1,1,2)

locs = [1,2,4,5,7,8,12,13,15,16,18,19];
means = R2_M(1:12,2);
stds = R2_STD(1:12,2);
plot_performance_bars(means, stds,locs,2,1,2)

% Now let's compare best lasso and RF performance in both RMSE and R2
locs = [1,2,3,4,6,7,8,9,11,12,13,14];
means = RMSE_M([7,8,13,14,9,10,15,16,11,12,17,18],2);
stds = RMSE_STD([7,8,13,14,9,10,15,16,11,12,17,18],2);
plot_performance_bars(means, stds,locs,1,2,2)

locs = [1,2,3,4,6,7,8,9,11,12,13,14];
means = R2_M([7,8,13,14,9,10,15,16,11,12,17,18],2);
stds = R2_STD([7,8,13,14,9,10,15,16,11,12,17,18],2);
plot_performance_bars(means, stds,locs,2,2,2)
