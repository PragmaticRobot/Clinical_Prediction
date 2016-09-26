function [ gg ] = loadmodres( ez )
%loadmodres loads model results coming from R using textscan
% output will be a cell array with each cell being the result of each
% cross-validation run, giving the prediction made in that run for each
% subject.
% each row is a subject
% each column is a cell array with a cross-validation run
fid = fopen(ez);
pit = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','delimiter',',');
fclose(fid);
gg = zeros(26,100);
for i = 1:100
    gg(:,i) = cell2mat(pit{i+1}(2:end));
end

end

