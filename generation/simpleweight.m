function [ weights ] = simpleweight( test, highlight )
%SIMPLEWEIGHT Summary of this function goes here
%   Detailed explanation goes here

weights = highlight.maxdist - vnorm(bsxfun(@minus,test,highlight.point));
weights = weights./sum(weights);

end

