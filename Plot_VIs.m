function [ ] = Plot_VIs( XX , titles)
%Plot_VIs Takes an 8-row vector and plots the feature importance
%   The input has to be an 8-row vector of feature ranks as integers
figure
clf

title('Rank stability across cross-validation repeats')

X_ax = 1:1:size(XX,2);
for i = 1:size(XX,1);
    subplot(4,2,i);
    plot(X_ax, XX(i,:), 'k-','LineWidth',1.5);
    ylabel(titles(i))
end

end

