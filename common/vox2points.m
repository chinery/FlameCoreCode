function [ points ] = vox2points( vox, mode, sample, denoise )
%VOX2POINTS Summary of this function goes here
%   Detailed explanation goes here
if nargin < 4
    denoise = 0;
end

if nargin < 3
    sample = 100;
end

if nargin < 2
    mode = 1;
end

if nargin < 1
    error('Need an input volume');
end

N = size(vox);
box = [1 N(2);
    1 N(1);
    1 N(3)];
[x,y,z] = meshgrid( box(1,1):box(1,2), box(2,1):box(2,2),box(3,1):box(3,2) );

in = vox > 0;

% if nargin == 3
%     sf = sample/sum(vox(:));
% else
%     sf = 100/max(vox(:));
% end
sf = sample;

if mode == 0
    points = [x(in) y(in) z(in)]';
elseif mode == 1 || mode == 2
    scaleddensities = ceil(vox * sf);
    
    scaleddensities(scaleddensities < denoise) = 0; % hack to denoise
    
    % some super smart person on stack overflow came up with this very speedy
    % implementation. I don't know how.
    % https://stackoverflow.com/questions/1975772/matlab-array-manipulation
    x = x(scaleddensities>0)';
    y = y(scaleddensities>0)';
    z = z(scaleddensities>0)';
    scaleddensities = scaleddensities(scaleddensities>0)';
    index = zeros(1,sum(scaleddensities));
    index([1 cumsum(scaleddensities(1:end-1))+1]) = 1;
    index = cumsum(index);
    points = [x(index); y(index); z(index)];
    
% %     in = scaleddensities > 0;% find data points
% %     points = zeros(3,sum(scaleddensities(:)));
% %     insert = 1;
% %     while any(in(:))
% %         newpoints = [x(in) y(in) z(in)]';
% %         a = insert;
% %         b = insert+size(newpoints,2)-1;
% %         points(:,a:b) = newpoints;
% %         insert = insert + size(newpoints,2) + 1;
% %         scaleddensities(scaleddensities > 0) = scaleddensities(scaleddensities > 0) - 1;
% %         in = scaleddensities > 0;
% %     end
    if mode == 2
        points = points + (rand(size(points)) -0.5);
    end
elseif mode == 3
    points = [x(in) y(in) z(in) vox(in)]';

    sig = std( points,[],2 ); % std of each dimension
    points = points ./ repmat( sig,1,size(points,2) ); % whiten
    points(1:4,:) = ((points(1:4,:)-min(min(points(1:4,:))))./(max(max(points(1:4,:)))-min(min(points(1:4,:)))))*100;
else
    error('Wrong mode');
end

points = bsxfun(@rdivide,points,size(vox)');


end

