% fixing coefficients so that they're reported on the standardized scale
% instead of the original scale

% coef * original = X * (original - mean)/sd
% X = coef*sd + mean/original
function [new_coefs] = fix_coefs(olds,N)
% olds are the "scaled" coefficients (Matrix of 100 x length(N))
% N is the vector of column numbers for the features
% load the whole table of "original"
Features = readtable('FeatureSet.csv');
Features(:,1) = [];

% Pick out features in N
feats = table2array(Features(:,N));
oldfeats = olds(:,N)