function [ points ] = cleanuppoints( points, vsize )
%CLEANUPPOINTS Summary of this function goes here
%   Detailed explanation goes here

[~,lab] = chincluster(points,vsize);
del = false(1,size(points,2));
for cl = 1:max(lab)
    part = points(:,lab==cl);
    if(size(part,2) < size(points,2)/100 || min(max(part,[],2)-min(part,[],2)) <= 4/vsize)
        del(lab==cl) = true;
    end
end
points(:,del) = []; 

end

