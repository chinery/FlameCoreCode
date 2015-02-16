function [clustCent,data2cluster,clusterVotes] = definedmeanshiftcluster(data,num,itr)
% Run mean shift clustering with various bandwidths until the right number of
% clusters is found.
% 
% data is data
% num is number of clusters to find

if(num == 1)
    data2cluster = ones(1,size(data,2));
    clustCent = mean(data,2);
    clusterVotes = ones(1,size(data,2));
    return;
end

if(nargin < 3)
    itr = 15;
end

d = size(data,1);

base = norm(max(data(1:d-1,:),[],2)-min(data(1:d-1,:),[],2));

bandwidths(1,1) = base/num;

i = 1;

while(true)
    b = bandwidths(i,1);
    [clustCent,data2cluster,clusterVotes] = CMeanShiftCluster(data,b);

    if(max(data2cluster) == num)
        return;
    elseif(i > itr)
        closeunder = max(bandwidths(bandwidths(:,2) < num,2));
        closeover = min(bandwidths(bandwidths(:,2) > num,2));
        error('Could not get this many clusters. Closest I got was %i or %i.',closeunder,closeover);
    end

    bandwidths(i,2) = max(data2cluster);
    
    bandwidths(i+1,1) = getnextband;
    i = i + 1;
end
    

    % nested functions
    function inval = invalid
        % if the latest result has a smaller bandwidth than a previous but a
        % smaller number of components (or larger/larger) then something is
        % terribly wrong.
        inval = any(bsxfun(@gt,bandwidths(i,1),bandwidths(1:i-1,1)) & ...
                    bsxfun(@gt,bandwidths(i,2),bandwidths(1:i-1,2))) || ...
                any(bsxfun(@lt,bandwidths(i,1),bandwidths(1:i-1,1)) & ...
                    bsxfun(@lt,bandwidths(i,2),bandwidths(1:i-1,2)));
    end
    
    function b = getnextband
        % if we've gone in the wrong direction
        if(i > 1 && ((bandwidths(i-1,2) < num && bandwidths(i,2) < bandwidths(i-1,2)) ||...
                    (bandwidths(i-1,2) > num && bandwidths(i,2) > bandwidths(i-1,2))))
            bandwidths(i-2:i-1,:) = [];
            i = i-2;
        elseif(invalid)
            % this is probably never going to work, but keep trying
            b = (bandwidths(i-1,1) + bandwidths(i,1))/2 + ...
                (rand*mean(bandwidths(:,1)) - rand*mean(bandwidths(:,1)))/10;
            return;
        end
        
        if(all(bandwidths(:,2) < num))
            b = base/(num+i-1);
        elseif(all(bandwidths(:,2) > num))
            if(num-i < 0)
                error('Error! num is likely implausibly low.');
            end
            b = base/(num-i+1);
        else
            % need to find the midpoint between the latest, closest results
            lt = find(bandwidths(:,2) < num);
            [~,ixlt] = min(num-bandwidths(lt,2));
            ixlt = find(bandwidths(:,2) == bandwidths(lt(ixlt),2),1,'last');

            gt = find(bandwidths(:,2) > num);
            [~,ixgt] = min(bandwidths(gt,2)-num);
            ixgt = find(bandwidths(:,2) == bandwidths(gt(ixgt),2),1,'last');

            b = (bandwidths(ixlt,1) + bandwidths(ixgt,1))/2;
        end
    end
end
