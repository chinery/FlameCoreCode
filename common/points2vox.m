function [ vox ] = points2vox( points, siz )
%POINTS2VOX Summary of this function goes here
%   Detailed explanation goes here

points = bsxfun(@times,points,siz');

points(1,:) = ceil(points(1,:));
points(2,:) = floor(points(2,:));
points(3,:) = ceil(points(3,:));
points(:,points(1,:)<1) = [];
points(:,points(2,:)<1) = [];
points(:,points(3,:)<1) = [];
points(:,points(1,:)>siz(1)) = [];
points(:,points(2,:)>siz(2)) = [];
points(:,points(3,:)>siz(3)) = [];


vox = zeros(siz);

for i = 1:length(points)
    vox(points(2,i), points(1,i), points(3,i)) = vox(points(2,i), points(1,i), points(3,i)) + 1;
end

% vox(:) = min(vox(:)./(mean(vox(vox(:)~=0))+10*std(vox(vox(:)~=0))),1);
vox(:) = vox(:)./max(vox(:));

end

