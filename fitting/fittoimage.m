clear all

twodflame.normm = [];
twodflame.norms = [];
twodflame.scale = [];
twodflame.unifo = [];
twodflame.maxim = [];
twodflame.maxes = [];

% a quick experiment...

img = double(imread('pinkflame.png'))./255;

img = imresize(img,2);

mask = rgb2gray(img)>0.045;

mask = imfilter(mask,fspecial('average',5));
% imshow(mask);

[x,y] = meshgrid(linspace(0,1,size(img,2)),linspace(1,0,size(img,1)));

grid = [x(:)';y(:)'];

points = grid(:,mask(:));

figure(3);
plot(points(1,:),points(2,:),'.');
axis equal

% [v, ~] = eig(cov(points'));
v = eye(2);

n = 40;
nn = (max(points(2,:))-min(points(2,:)))./(n-1);
qpoints = v'*points;
qpoints(2,:) = round(points(2,:)./nn).*nn;

[a,~] = count_unique(qpoints(2,:));
core = zeros(2,length(a));
core(2,:) = a';
core(1,:) = arrayfun(@(x)mean(qpoints(1,qpoints(2,:)==x)),unique(qpoints(2,:)));

core = v*core;

x = core(1,:);
y = core(2,:);

t = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2)])';
len = max(t);
t = t./len;

t = t.*100;

xfit = fit(t',x','smoothingspline','SmoothingParam',0.001);
yfit = fit(t',y','smoothingspline','SmoothingParam',0.001);

bestcore = [];
bestcore(1,:) = xfit(t)';
bestcore(2,:) = yfit(t)';

bestcore(:,end) = [];

hold on;
plot(bestcore(1,:),bestcore(2,:),'xk');

%% fidn parameters
r = img(:,:,1);
g = img(:,:,2);
b = img(:,:,3);
grey = rgb2gray(img);
% density.r = r(mask(:))';
% density.g = g(mask(:))';
% density.b = b(mask(:))';
% density.grey = grey(mask(:))';
mask = imfilter(double(mask),fspecial('average',10))>0;
points = grid(:,mask(:));
density.r = r(mask(:))';
density.g = g(mask(:))';
density.b = b(mask(:))';
density.grey = grey(mask(:))';


ival = 40;
bestcore = interp1(bestcore',linspace(1,size(bestcore,2),ival),'spline','extrap')';

threedcore(1,:) = bestcore(1,:);
threedcore(2,:) = 0.5;
threedcore(3,:) = bestcore(2,:);
twodflame.cores{1} = threedcore;

[dist2core, point2core] = pdist2(bestcore',points','euclidean','smallest',1);
vec2core = vnorm(points - bestcore(:,point2core));

[a,b] = count_unique(point2core);
range = a;

if(length(range) ~= length(bestcore))
    bestcore = bestcore(:,range);
    bestcore = interp1(bestcore',linspace(1,size(bestcore,2),ival),'spline','extrap')';
    
    [dist2core, point2core] = pdist2(bestcore',points','euclidean','smallest',1);
    vec2core = vnorm(points - bestcore(:,point2core));
end

[a,b] = count_unique(point2core);
range = a;

paramfit = arrayfun(@(x)twodestimateparameters(dist2core(point2core==x)',density.grey(point2core==x)),range);

a = range;
normm = zeros(1,length(bestcore));
norms = zeros(1,length(bestcore));
unifo = zeros(1,length(bestcore));
scale = zeros(1,length(bestcore));
maxim = zeros(1,length(bestcore));

normm(a) = arrayfun(@(x)x.mu,paramfit);
norms(a) = arrayfun(@(x)x.sig,paramfit);
unifo(a) = arrayfun(@(x)x.u,paramfit);
scale(a) = arrayfun(@(x)x.s,paramfit);
maxim(a) = arrayfun(@(x)x.m,paramfit);

twodflame.normm{1}(1:3,:) = repmat(normm,3,1);
twodflame.norms{1}(1:3,:) = repmat(norms,3,1);
twodflame.maxim{1}(1:3,:) = repmat(maxim,3,1);

for col = 3:-1:1
    if(col == 1)
        dens = density.r;
    elseif(col == 2)
        dens = density.g;
    else
        dens = density.b;
    end
        
    
    scalefit = arrayfun(@(x)twodfitscaleparameter(dist2core(point2core==x)',dens(point2core==x),...
                normm(x),norms(x),unifo(x),scale(x),maxim(x)),range);
            
    a = range;
    cunifo = zeros(1,length(bestcore));
    cscale = zeros(1,length(bestcore));
    cunifo(a) = arrayfun(@(x)x.u,scalefit);
    cscale(a) = arrayfun(@(x)x.s,scalefit);

    twodflame.unifo{1}(col,:) = cunifo;
    twodflame.scale{1}(col,:) = cscale;
    
%     maxes = zeros(length(bestcore),1);
%     maxes(a) = arrayfun(@(x)max(dist2core(point2core==x)),range);
% 
%     twodflame.maxes{1}(col,:) = maxes;
    
end

% get rid of the bottom
rsize = @(x)(interp1(x',linspace(2,size(x,2),ival))');
twodflame.cores{1} = rsize(twodflame.cores{1});
twodflame.normm{1} = rsize(twodflame.normm{1});
twodflame.norms{1} = rsize(twodflame.norms{1});
twodflame.scale{1} = rsize(twodflame.scale{1});
twodflame.unifo{1} = rsize(twodflame.unifo{1});
twodflame.maxim{1} = rsize(twodflame.maxim{1});



figure;displaysingleflame(twodflame);