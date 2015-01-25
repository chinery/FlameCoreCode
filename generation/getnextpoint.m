function [newpoint,direction] = getnextpoint(currentpoint,direction,trainingdata,matching,sigma,distscale,densityfn,weightfn,parameters)

[~, map] = ismember(matching(1,:),1:size(trainingdata,2));

[dist,ix] = pdist2(trainingdata(:,matching(1,:))',currentpoint','euclidean','smallest',parameters.nearest.num);
ix = map(ix);
dist = dist./distscale;

if(parameters.nearest.num > 1 && any(dist < parameters.nearest.distance) && any(dist > parameters.nearest.distance))
    ix(dist > parameters.nearest.distance) = [];
    dist(dist > parameters.nearest.distance) = [];
end

dist = dist(1:min(parameters.nearest.num,length(dist)));
ix = ix(1:min(parameters.nearest.num,length(ix)));

nextix = matching(2,ismember(matching(1,:),ix));
nextpoints = trainingdata(:,nextix);

if(length(ix) > 1)
    if(length(currentpoint) > 1 && any(direction ~= 0))
        potentialdirections = bsxfun(@minus,nextpoints,currentpoint);
        angles = vectorangle(potentialdirections,repmat(direction,1,size(potentialdirections,2)));
        weights = 1-(angles./pi);
    else
        dist(dist == 0) = parameters.nearest.exactdistcorrection;
        weights = 1./dist;
    end
    weights = weights(:)./sum(weights(:));
    if(~isempty(weightfn))
        weights = weights'.*weightfn(nextpoints);
        weights = weights./sum(weights);
    end
    point = randsample(1:length(ix),1,true,weights);
else
    point = 1;
end

nextpoint = nextpoints(:,point);

% % % This code will move in the same direction from the current point
% % % rather than snapping back to the training data
% % % Note: untested
% nearpoint = matchingdata(:,ix(point));
% vec = nextpoint-nearpoint;
% nextpoint = currentpoint + vec;


newpoints = mvnrnd(nextpoint',parameters.gaussian.scale*sigma',parameters.gaussian.num)';
if(parameters.gaussian.scale == 0)
    likelihood = ones(parameters.gaussian.num,1);
else
    likelihood = mvnpdf(newpoints',nextpoint',parameters.gaussian.scale*sigma')';
end

if(any(direction ~= 0))
    potentialdirections = bsxfun(@minus,newpoints,currentpoint);
    angles = vectorangle(potentialdirections,repmat(direction,1,size(potentialdirections,2)));
    angleweights = 1-(angles./pi);
    angleweights(isnan(angleweights)) = 0.5;
else
    angleweights = ones(1,parameters.gaussian.num);
end
angleweights = angleweights./sum(angleweights);

likelihood = likelihood./sum(likelihood);
dest = densityfn(newpoints);
dest = dest./sum(dest);

weights = angleweights .* likelihood .* dest;
weights = weights./sum(weights);

newpoint = newpoints(:,randsample(1:size(newpoints,2),1,true,weights));

direction = newpoint - currentpoint;

end