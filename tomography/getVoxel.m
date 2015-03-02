function [ xmin, xmax, ymin, ymax, zmin, zmax ] = getVoxel( number, totalVoxelsPerDimension )
%GETVOXEL returns boundaries for given voxel number
%   voxels are numbered in order of x, then y, then z. So if you ask for
%   voxel 1 you get the one starting at the gobal voxelfrom, if you ask for
%   2 you get the one that is one along in the x direction.

%voxels contained between
global voxelfrom 
global voxelto
a = voxelfrom;
b = voxelto;

% in tests using ind2sub was significantly slower than doing it manually.
znum = floor((number-1)/totalVoxelsPerDimension^2)+1;
remainder = number - (znum-1)*totalVoxelsPerDimension^2;
ynum = floor((remainder-1)/totalVoxelsPerDimension)+1;
xnum = remainder - (ynum-1)*totalVoxelsPerDimension;

interval = (b-a)./totalVoxelsPerDimension;
offset = a;
xmax = xnum * interval(1) + offset(1);
ymax = ynum * interval(2) + offset(2);
zmax = znum * interval(3) + offset(3);

xmin = xmax - interval(1);
ymin = ymax - interval(2);
zmin = zmax - interval(3);


end


