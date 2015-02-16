function [ smoothed ] = surfacesmooth( u1, u2, data, window, overlap, weight, include )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if(nargin < 6 || isempty(weight))
    weight = ones(1,size(data,2));
end

if(nargin < 7)
    include = false;
end

weight = repmat(weight,size(data,1),1);

if(include)
    exclude = false(size(data));
    exclude(weight==0) = true;
else
    lasts = false(size(data));
    lasts(end,:) = true;
    exclude = data == 0 | lasts;
end

weight(weight==0) = 0.1;

offset = window-overlap;
startidx = 1:offset:size(data,2)-window+1;
endidx = startidx + window - 1;

if(endidx(end) ~= size(data,2))
    startidx(end+1) = size(data,2)-window+1;
    endidx(end+1) = size(data,2);
end

% startidx = 1:offset:100-window+1;
% endidx = startidx + window - 1;

agg = zeros(size(data));
tot = zeros(size(data));

progressbar(0);

for i = 1:length(startidx)
    from = startidx(i);
    to = endidx(i);
    
    t1 = u1(:,from:to);
    t2 = u2(:,from:to);
    d = data(:,from:to);
    ex = exclude(:,from:to);
    we = weight(:,from:to);
    
    if(all(ex))
        continue;
    end
    
    dfit = fit([t1(:) t2(:)],d(:),'lowess','Exclude',ex(:),'Weights',we(:));%,'Normalize','on');

    dfitsamp = dfit(t1(:),t2(:));

    dfitsamp = reshape(dfitsamp,size(data,1),[]);
    
    agg(:,from:to) = agg(:,from:to) + 1;
    tot(:,from:to) = tot(:,from:to) + dfitsamp;
    
    progressbar(i/length(startidx));
end

smoothed = tot./agg;

smoothed(exclude) = data(exclude);

% result = indata;
% result(2:end-1,:) = smoothed;
% 
% smoothed = result;

end

