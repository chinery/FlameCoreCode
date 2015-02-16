function [ centres, label, votes ] = chincluster( points, maxvoxdim )
%CHINCLUSTER Summary of this function goes here
%   Detailed explanation goes here

data = unique(points','rows')';

smalllabel = zeros(1,length(data));
currentcluster = 1;

while(any(smalllabel == 0))
    unassigned = find(smalllabel == 0);
    left = data(:,smalllabel == 0);
    
    ix = randi(size(left,2));
    
    dists = vnorm(bsxfun(@minus,left,left(:,ix)));

    ndist = 1.5/maxvoxdim;

    neighbours = find(dists < ndist);

    if(isempty(neighbours))
        smalllabel(unassigned(ix)) = currentcluster;
        currentcluster = currentcluster+1;
        continue;
    end
    
    totalix = ix;
    while(true)
        totalix = unique([totalix neighbours]);
        ix = neighbours;

        rest = left;
        rest(:,totalix) = [];
        restix = 1:size(left,2);
        restix(totalix) = [];

        [D,I] = pdist2(rest', left(:,ix)','euclidean','Smallest',19);

        neighbours = unique(restix(I(D < ndist)));
        
        if(isempty(neighbours))
            smalllabel(unassigned(totalix)) = currentcluster;
            currentcluster = currentcluster+1;
            break;
        end
    end
end

[~,index] = ismember(points',data','rows');
label = smalllabel(index);

clusters = currentcluster-1;
centres = zeros(3,clusters);
votes = zeros(clusters,size(points,2));
for i = 1:clusters
    centres(:,i) = mean(points(:,label==i),2);
    votes(i,:) = label==i;
end


end

