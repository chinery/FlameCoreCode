function [ num ] = guesslicks( points, vsize, cut )
%GUESSLICKS Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 3)
    cut = true;
end

if(cut)
    points(:,points(3,:) < 0.2) = [];
end

[~,lab] = chincluster(points,vsize);

valid = 0;
for cl = 1:max(lab)
    part = points(:,lab==cl);
    if(~(size(part,2) < size(points,2)/100 || min(max(part,[],2)-min(part,[],2)) <= 4/vsize))
        valid = valid+1;
    end
end

num = valid;

end

