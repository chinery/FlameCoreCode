function [core, flat] = fitcore(points,lrgpoints,maxvoxdim,oldmethod,rotate)

%
% points    the points to fit to
% lrgpoints gives points that are used in the averaging process, but not used to
%           generate new points. E.g., if one lick of flame is connected to a
%           base, then input the flame + base points as lrgpoints. This will
%           give a smooth connection between the lick and the base, even though
%           no 'flame points' will lie within the base.
%           If you leave it blank, it will create a block of points around the
%           base of the flame since that is the most problimatic area (don't
%           want the core to drive sharply into the corner, but to stay centred)
% maxvoxdim the resolution of the reconstruction, e.g. 128 for 128x128x128
% oldmethod whether to use the older sliced based method that is more robust to
%           noise, default false, but you should compare both for your data
% rotate    whether to rotate to align with principle directions, default true
%           (only applicable if oldmethod = true)



if(nargin<5)
    rotate = true;
end

global nn
nn = 1/maxvoxdim;
doplot = false;

% parameters
num = 1.5;
gap = 1;

% for old time's sake
[v, ~] = eig(cov(points'));
for k = size(v,2):-1:1
    angle(k) = atan2(norm(cross(v(:,k),[0 0 1])),dot(v(:,k),[0 0 1]));
end
[~, k] = min(angle);

cen = mean(points,2);
flat = k ~= 3 & cen(3) <= 20/maxvoxdim;

if(isempty(lrgpoints))
    upmost = v(:,k);
    v(:,k) = [];
    v = [v upmost];

    rpoints = v'*points;
    
    [x,y,z] = meshgrid(min(rpoints(1,:)):nn:max(rpoints(1,:)),min(rpoints(2,:)):nn:max(rpoints(2,:)),(min(rpoints(3,:))-5*nn):nn:(min(rpoints(3,:)+40*nn)));
    qlrpoints = [rpoints repmat([x(:)';y(:)';z(:)'],1,10)];%part;%points(:,ix==i | ix==max(ix));
    lrgpoints = v*qlrpoints;
end

% prune
lrgpoints(:,lrgpoints(1,:) < min(points(1,:))) = [];
lrgpoints(:,lrgpoints(2,:) < min(points(2,:))) = [];
lrgpoints(:,lrgpoints(1,:) > max(points(1,:))) = [];
lrgpoints(:,lrgpoints(2,:) > max(points(2,:))) = [];



if(exist('oldmethod','var') && oldmethod)
    core = rollbackcore(points,lrgpoints,maxvoxdim,rotate);
    if(~isempty(core) && size(core,2) > 2)
        core = splinesmooth(core,points,4);
    else
        core = [];
    end
    return
end

% proper stuff

points = unique(points','rows')';
lrgpoints = unique(lrgpoints','rows')';

points(3,:) = round(points(3,:)./nn)*nn;

middlez = (max(points(3,:)) - min(points(3,:)))/2 + min(points(3,:));

disc = points(:,points(3,:) <= middlez+num*nn & points(3,:) >= middlez-num*nn);
numel(unique(disc(3,:)));
middle = mean(disc,2);

centop = mean(disc(:,disc(3,:) == max(disc(3,:))),2);
cenbot = mean(disc(:,disc(3,:) == min(disc(3,:))),2);

upwardsness = norm(centop(1:2) - cenbot(1:2));
%
% figure; hold on;
% plot3(points(1,:),points(2,:),points(3,:),'.','markersize',1);
% plot3(disc(1,:),disc(2,:),disc(3,:),'xr');
% plot3(middle(1),middle(2),middle(3),'*g');
% plot3(centop(1),centop(2),centop(3),'*k');
% plot3(cenbot(1),cenbot(2),cenbot(3),'*k');

R = eye(3);
while(true)
    u = centop - cenbot;
    v = [0 0 1]';
    
    u = u./norm(u);
    tmpR = vecRotMat(u',v');
    
    rotpoints = tmpR*R*points;
    rotpoints(3,:) = round(rotpoints(3,:)./nn)*nn;
    
    rmiddle = tmpR*R*middle;
    rmiddle(3) = round(rmiddle(3)./nn)*nn;
    rmiddlez = rmiddle(3);
    
    rdisc = rotpoints(:,rotpoints(3,:) <= rmiddlez+num*nn & rotpoints(3,:) >= rmiddlez-num*nn);
    numel(unique(rdisc(3,:)));
    
    try
        [~,inx] = pruneplane(rdisc(1:2,:),rmiddle(1:2),maxvoxdim);
        rdisc = rdisc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    centop = mean(rdisc(:,rdisc(3,:) == max(rdisc(3,:))),2);
    cenbot = mean(rdisc(:,rdisc(3,:) == min(rdisc(3,:))),2);
    
    newup = norm(centop(1:2) - cenbot(1:2));
    
    if(newup < upwardsness)
        R = tmpR*R;
        upwardsness = newup;
        
        test = R*[0 0 1]';
        angle = atan2(norm(cross(test,[0 0 1])),dot(test,[0 0 1]));
        if(angle/pi > 0.3)
            R = eye(3);
            warning('I don''t like where this is going'); %#ok
            break;
        end
        
        %         figure; hold on;
        %         plot3(rotpoints(1,:),rotpoints(2,:),rotpoints(3,:),'.','markersize',1);
        %         plot3(rdisc(1,:),rdisc(2,:),rdisc(3,:),'xr');
        %         plot3(rmiddle(1),rmiddle(2),rmiddle(3),'*g');
        %
        %         plot3(centop(1),centop(2),centop(3),'*k');
        %         plot3(cenbot(1),cenbot(2),cenbot(3),'*k');
    else
        break;
    end
end


rotpoints = R*points;
rotpoints(3,:) = round(rotpoints(3,:)./nn)*nn;

rmiddle = R*middle;
rmiddle(3) = round(rmiddle(3)./nn)*nn;
rmiddlez = rmiddle(3);

% backup the current state of things for when we do the lower half
halfwayR = R;
originalmiddle = middle;

% find the first point
rdisc = rotpoints(:,rotpoints(3,:) <= rmiddlez+num*nn & rotpoints(3,:) >= rmiddlez-num*nn);

try
    [~,inx] = pruneplane(rdisc(1:2,:),rmiddle(1:2),maxvoxdim);
    rdisc = rdisc(:,inx);
catch MExc
    if(~strcmp(MExc.identifier(1:5),'Prune'))
        rethrow(MExc);
    end
end

rmiddle = mean(rdisc,2);

top(:,1) = R'*rmiddle;

lrgrpoints = R*lrgpoints;
lrgrpoints(3,:) = round(lrgrpoints(3,:)./nn)*nn;

%
abort = false;
babort = false;
oob = false;

%% do the top half
c = 2;
baseR = R;  %baseR goes from points to the halfway rotated points
while(true)
    %move up the flame
    middlez = rmiddle(3) + gap*nn + num*nn;
    disc = rotpoints(:,rotpoints(3,:) <= middlez+num*nn & rotpoints(3,:) >= middlez-num*nn);
    
    % when going up the flame, stop before reaching the very tip to get the
    % effect of a 'ball' at the top
    D = pdist2(mean(disc,2)',disc','euclidean','Largest',10);
    [~,ix] = max(rotpoints(3,:));
    if(mean(D) > norm(rotpoints(:,ix) - mean(disc,2)))
        break;
    end
    
    try
        [~,inx] = pruneplane(disc(1:2,:),rmiddle(1:2),maxvoxdim);
        disc = disc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    if(size(disc,2) < 5 || numel(unique(disc(3,:))) == 1)
        break;
    end
    
    disc = lrgrpoints(:,lrgrpoints(3,:) <= middlez+num*nn & lrgrpoints(3,:) >= middlez-num*nn);
    
    try
        [~,inx] = pruneplane(disc(1:2,:),rmiddle(1:2),maxvoxdim);
        disc = disc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    middle = mean(disc,2);
    
    centop = mean(disc(:,disc(3,:) == max(disc(3,:))),2);
    cenbot = mean(disc(:,disc(3,:) == min(disc(3,:))),2);
    
    upwardsness = norm(centop(1:2) - cenbot(1:2));
    
    if(doplot)
        figure(10);clf; hold on;
        plot3(rotpoints(1,:),rotpoints(2,:),rotpoints(3,:),'.','markersize',1);
        plot3(disc(1,:),disc(2,:),disc(3,:),'xr');
        plot3(middle(1),middle(2),middle(3),'*g');
        plot3(centop(1),centop(2),centop(3),'*k');
        plot3(cenbot(1),cenbot(2),cenbot(3),'*k');
        view(20,15);
        pause(0.2)
%         %         blackplot
%         blahblah = blahblah + 1;
%         print('-dpng',sprintf('./results/fitting2/%.4i',blahblah));
    end
    
    newR = eye(3); %this is this modification we're looking to make
    while(true)
        u = centop - cenbot;
        v = [0 0 1]';
        
        u = u./norm(u);
        tmpR = vecRotMat(u',v');
        
        rotpoints = tmpR*newR*baseR*points;
        rotpoints(3,:) = round(rotpoints(3,:)./nn)*nn;
        
        lrgrpoints = tmpR*newR*baseR*lrgpoints;
        lrgrpoints(3,:) = round(lrgrpoints(3,:)/nn)*nn;
        
        rmiddle = tmpR*newR*middle;
        rmiddle(3) = round(rmiddle(3)./nn)*nn;
        rmiddlez = rmiddle(3);
        
        rdisc = lrgrpoints(:,lrgrpoints(3,:) <= rmiddlez+num*nn & lrgrpoints(3,:) >= rmiddlez-num*nn);
        
        try
            [~,inx] = pruneplane(rdisc(1:2,:),rmiddle(1:2),maxvoxdim);
            rdisc = rdisc(:,inx);
        catch MExc
            if(strcmp(MExc.identifier(1:5),'Prune'))
                break;
            else
                rethrow(MExc);
            end
        end
        
        centop = mean(rdisc(:,rdisc(3,:) == max(rdisc(3,:))),2);
        cenbot = mean(rdisc(:,rdisc(3,:) == min(rdisc(3,:))),2);
        
        indisc = rotpoints(:,rotpoints(3,:) <= rmiddlez+num*nn & rotpoints(3,:) >= rmiddlez-num*nn);
        
        try
            [~,inx] = pruneplane(indisc(1:2,:),rmiddle(1:2),maxvoxdim);
            indisc = indisc(:,inx);
        catch MExc
            if(strcmp(MExc.identifier(1:5),'Prune'))
                break;
            else
                rethrow(MExc);
            end
        end
        
        if(size(indisc,2) < 5 || numel(unique(indisc(3,:))) == 1)
            newR = eye(3);
            break;
        end
        
        
        if(doplot)
            figure(20);clf; hold on;
            plot3(rotpoints(1,:),rotpoints(2,:),rotpoints(3,:),'.','markersize',1);
            plot3(rdisc(1,:),rdisc(2,:),rdisc(3,:),'xr');
            plot3(rmiddle(1),rmiddle(2),rmiddle(3),'*g');
            plot3(centop(1),centop(2),centop(3),'*k');
            plot3(cenbot(1),cenbot(2),cenbot(3),'*k');
            view(20,15);
            pause(0.2)
        end
        
        
        newup = norm(centop(1:2) - cenbot(1:2));
        
        if(newup < upwardsness)
            newR = tmpR*newR;
            upwardsness = newup;
            
            test = newR*[0 0 1]';
            angle = atan2(norm(cross(test,[0 0 1])),dot(test,[0 0 1]));
            if(angle/pi > 0.3)
                break;
            end
        else
            break;
        end
    end
    
    rotpoints = newR*baseR*points;
    rotpoints(3,:) = round(rotpoints(3,:)./nn)*nn;
    
    lrgrpoints = newR*baseR*lrgpoints;
    lrgrpoints(3,:) = round(lrgrpoints(3,:)/nn)*nn;
    
    rmiddle = newR*middle;
    rmiddle(3) = round(rmiddle(3)./nn)*nn;
    rmiddlez = rmiddle(3);
    
    rdisc = lrgrpoints(:,lrgrpoints(3,:) <= rmiddlez+num*nn & lrgrpoints(3,:) >= rmiddlez-num*nn);
    
    try
        [~,inx] = pruneplane(rdisc(1:2,:),rmiddle(1:2),maxvoxdim);
        rdisc = rdisc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    rmiddle = mean(rdisc,2);
    
    indisc = rotpoints(:,rotpoints(3,:) <= rmiddlez+num*nn & rotpoints(3,:) >= rmiddlez-num*nn);
    
    try
        [~,inx] = pruneplane(indisc(1:2,:),rmiddle(1:2),maxvoxdim);
        indisc = indisc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    if(size(indisc,2) < 5)
        break
    end
    
    bbox = [min(indisc,[],2) max(indisc,[],2)];
    
    if(any(rmiddle(1:2) < bbox(1:2,1)) || any(rmiddle(1:2) > bbox(1:2,2)))
        break;
    end
    
    baseR = newR*baseR;
    
    top(:,c) = baseR'*rmiddle;
    
    rtop = baseR*top;
    if(any(rtop(3,c) < rtop(3,:) - nn))
        % if this latest point is lower than any of the previous
        %         if(~any(vnorm(bsxfun(@minus,top(1:2,top(3,c) < top(:,c)), top(1:2,c))) <...
        %                 vnorm(bsxfun(@minus,top(3,top(3,c) < top(:,c)), top(3,c)))))
        %             % and not 'next door'
        %             abort = true;
        %         end
        abort = true;
        break
    end
    
    c = c+1;
end

if(doplot)
    figure(10); hold on;
    plot3(points(1,:),points(2,:),points(3,:),'.','markersize',1);
    plot3(top(1,:),top(2,:),top(3,:),'k*');
    % blackplot
%     blahblah = blahblah + 1;
%     print('-dpng',sprintf('./results/fitting2/%.4i',blahblah));
end

% something in here to check the fit of the top half?

% now do the bottom half
R = halfwayR;
rotpoints = R*points;
rotpoints(3,:) = round(rotpoints(3,:)./nn)*nn;

rmiddle = R*originalmiddle;
rmiddle(3) = round(rmiddle(3)./nn)*nn;
rmiddlez = rmiddle(3);

rdisc = rotpoints(:,rotpoints(3,:) <= rmiddlez+num*nn & rotpoints(3,:) >= rmiddlez-num*nn);

try
    [~,inx] = pruneplane(rdisc(1:2,:),rmiddle(1:2),maxvoxdim);
    rdisc = rdisc(:,inx);
catch MExc
    if(~strcmp(MExc.identifier(1:5),'Prune'))
        rethrow(MExc);
    end
end

rmiddle = mean(rdisc,2);

lrgrpoints = R*lrgpoints;
lrgrpoints(3,:) = round(lrgrpoints(3,:)./nn)*nn;

bottom = [];
c = 1;
baseR = R;  %baseR goes from points to the halfway rotated points
while(true)
    %     if(abort)
    %         break;
    %     end
    
    %move down the flame
    middlez = rmiddle(3) - gap*nn - num*nn;
    disc = rotpoints(:,rotpoints(3,:) <= middlez+num*nn & rotpoints(3,:) >= middlez-num*nn);
    
    try
        [~,inx] = pruneplane(disc(1:2,:),rmiddle(1:2),maxvoxdim);
        disc = disc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    if(size(disc,2) < 5 || numel(unique(disc(3,:))) == 1)
        break;
    end
    
    disc = lrgrpoints(:,lrgrpoints(3,:) <= middlez+num*nn & lrgrpoints(3,:) >= middlez-num*nn);
    
    try
        [~,inx] = pruneplane(disc(1:2,:),rmiddle(1:2),maxvoxdim);
        disc = disc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    middle = mean(disc,2);
    
    centop = mean(disc(:,disc(3,:) == max(disc(3,:))),2);
    cenbot = mean(disc(:,disc(3,:) == min(disc(3,:))),2);
    
    upwardsness = norm(centop(1:2) - cenbot(1:2));
    
    if(doplot)
        figure(10);clf; hold on;
        plot3(rotpoints(1,:),rotpoints(2,:),rotpoints(3,:),'.','markersize',1);
        plot3(disc(1,:),disc(2,:),disc(3,:),'xr');
        plot3(middle(1),middle(2),middle(3),'*g');
        plot3(centop(1),centop(2),centop(3),'*k');
        plot3(cenbot(1),cenbot(2),cenbot(3),'*k');
        view(20,15);
        %     blackplot
%         blahblah = blahblah + 1;
%         print('-dpng',sprintf('./results/fitting/%.4i',blahblah));
    end
    
    newR = eye(3); %this is this modification we're looking to make
    while(true)
        u = centop - cenbot;
        v = [0 0 1]';
        
        %             u1 = u/norm(u);
        %             v1 = v/norm(v);
        %
        %             k = cross(u1,v1);
        %             % Rodrigues's formula:
        %             costheta = dot(u1,v1);
        %             newR =[ 0 -k(3) k(2);
        %                  k(3) 0 -k(1);
        %                 -k(2) k(1) 0];
        %             newR = costheta*eye(3) + newR + k*k'*(1-costheta)/sum(k.^2);
        
        u = u./norm(u);
        tmpR = vecRotMat(u',v');
        
        rotpoints = tmpR*newR*baseR*points;
        rotpoints(3,:) = round(rotpoints(3,:)./nn)*nn;
        
        lrgrpoints = tmpR*newR*baseR*lrgpoints;
        lrgrpoints(3,:) = round(lrgrpoints(3,:)/nn)*nn;
        
        rmiddle = tmpR*newR*middle;
        rmiddle(3) = round(rmiddle(3)./nn)*nn;
        rmiddlez = rmiddle(3);
        
        rdisc = lrgrpoints(:,lrgrpoints(3,:) <= rmiddlez+num*nn & lrgrpoints(3,:) >= rmiddlez-num*nn);
        
        try
            [~,inx] = pruneplane(rdisc(1:2,:),rmiddle(1:2),maxvoxdim);
            rdisc = rdisc(:,inx);
        catch MExc
            if(strcmp(MExc.identifier(1:5),'Prune'))
                break;
            else
                rethrow(MExc);
            end
        end
        
        centop = mean(rdisc(:,rdisc(3,:) == max(rdisc(3,:))),2);
        cenbot = mean(rdisc(:,rdisc(3,:) == min(rdisc(3,:))),2);
        
        indisc = rotpoints(:,rotpoints(3,:) <= rmiddlez+num*nn & rotpoints(3,:) >= rmiddlez-num*nn);
        
        try
            [~,inx] = pruneplane(indisc(1:2,:),rmiddle(1:2),maxvoxdim);
            indisc = indisc(:,inx);
        catch MExc
            if(strcmp(MExc.identifier(1:5),'Prune'))
                break;
            else
                rethrow(MExc);
            end
        end
        
        
        if(size(indisc,2) < 5 || numel(unique(indisc(3,:))) == 1)
            newR = eye(3);
            break;
        end
        
        if(doplot)
            figure(20);clf; hold on;
            plot3(rotpoints(1,:),rotpoints(2,:),rotpoints(3,:),'.','markersize',1);
            plot3(rdisc(1,:),rdisc(2,:),rdisc(3,:),'xr');
            plot3(rmiddle(1),rmiddle(2),rmiddle(3),'*g');
            plot3(centop(1),centop(2),centop(3),'*k');
            plot3(cenbot(1),cenbot(2),cenbot(3),'*k');
            view(20,15);
        end
        
        newup = norm(centop(1:2) - cenbot(1:2));
        
        if(newup < upwardsness)
            newR = tmpR*newR;
            upwardsness = newup;
            
            test = newR*[0 0 1]';
            angle = atan2(norm(cross(test,[0 0 1])),dot(test,[0 0 1]));
            if(angle/pi > 0.3)
                %                 newR = eye(3);
                %                 warning('I don''t like where this is going'); %#ok
                break;
            end
        else
            break;
        end
    end
    
    rotpoints = newR*baseR*points;
    rotpoints(3,:) = round(rotpoints(3,:)./nn)*nn;
    
    lrgrpoints = newR*baseR*lrgpoints;
    lrgrpoints(3,:) = round(lrgrpoints(3,:)/nn)*nn;
    
    rmiddle = newR*middle;
    rmiddle(3) = round(rmiddle(3)./nn)*nn;
    rmiddlez = rmiddle(3);
    
    rdisc = lrgrpoints(:,lrgrpoints(3,:) <= rmiddlez+num*nn & lrgrpoints(3,:) >= rmiddlez-num*nn);
    
    try
        [~,inx] = pruneplane(rdisc(1:2,:),rmiddle(1:2),maxvoxdim);
        rdisc = rdisc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    rmiddle = mean(rdisc,2);
    
    indisc = rotpoints(:,rotpoints(3,:) <= rmiddlez+num*nn & rotpoints(3,:) >= rmiddlez-num*nn);
    
    try
        [~,inx] = pruneplane(indisc(1:2,:),rmiddle(1:2),maxvoxdim);
        indisc = indisc(:,inx);
    catch MExc
        if(strcmp(MExc.identifier(1:5),'Prune'))
            break;
        else
            rethrow(MExc);
        end
    end
    
    if(size(indisc,2) < 5)
        break;
    end
    bbox = [min(indisc,[],2) max(indisc,[],2)];
    
    if(any(rmiddle(1:2) < bbox(1:2,1)) || any(rmiddle(1:2) > bbox(1:2,2)))
        break;
    end
    
    baseR = newR*baseR;
    
    bottom(:,c) = baseR'*rmiddle;
    
    rbottom = baseR*bottom;
    rtop = baseR*top;
    if(any(rbottom(3,c) > rbottom(3,:) + nn) || any(rbottom(3,c) > rtop(3,:) + nn))
        % if this latest point is higher than any of the previous
        %         if(~any(vnorm(bsxfun(@minus,bottom(1:2,bottom(3,c) < bottom(:,c)), bottom(1:2,c))) <...
        %                 vnorm(bsxfun(@minus,bottom(3,bottom(3,c) < bottom(:,c)), bottom(3,c)))))
        %             % and not 'next door'
        %             abort = true;
        %         end
        babort = true;
        break
    end
    c = c+1;
end

if(doplot)
    figure(10)
    plot3(bottom(1,:),bottom(2,:),bottom(3,:),'k*');
    % blackplot
%     blahblah = blahblah + 1;
%     print('-dpng',sprintf('./results/fitting/%.4i',blahblah));
end

core = [fliplr(bottom) top];

[regcore, err] = splinesmooth(core,points,4);
samp = ppval(regcore,linspace(0,1,100));
bbox = [min(points(1:2,:),[],2) max(points(1:2,:),[],2)];
if(sum(sum(bsxfun(@lt,samp(1:2,:),bbox(:,1)) | bsxfun(@gt,samp(1:2,:),bbox(:,2)))) > 10)
    oob = true;
    abort = true;
end

if(size(core,2) < 3 || abort || babort)
    
    if(size(core,2) < 3)
        core = rollbackcore(points,lrgpoints,maxvoxdim);
        [core1, err1] = splinesmooth(core,points,4);
        core = rollbackcore(points,points,maxvoxdim);
        [core2, err2] = splinesmooth(core,points,4);
        if(err2 < err1)
            core = core2;
            err = err2;
        else
            core = core1;
            err = err1;
        end
        warning('rollback')
    elseif(abort || babort)
        [testcore,testerr] = splinesmooth(core,points,4);
        
        if(abort)
            core = rollbackcore(points,lrgpoints,maxvoxdim);
            [core1, err1] = splinesmooth(core,points,4);
            core = rollbackcore(points,points,maxvoxdim);
            [core2, err2] = splinesmooth(core,points,4);
            if(err2 < err1)
                core = core2;
                err = err2;
            else
                core = core1;
                err = err1;
            end
        else
            botpoints = points(:,points(3,:) < originalmiddle(3));
            lbotpoints = lrgpoints(:,lrgpoints(3,:) < originalmiddle(3));
            
            core = [rollbackcore(botpoints,lbotpoints,maxvoxdim) top];
            [core, err] = splinesmooth(core,points,4);
        end
        
        if(testerr < err && ~oob)
            core = testcore;
            err = testerr;
        else
            warning('rollback')
        end
        
    end
else
    [core,err] = splinesmooth(core,points,4);
end

if(doplot)
    num = ((max(points(3,:))-min(points(3,:)))./nn) + 1;
    samp = ppval(core,linspace(0,1,num));
    figure; hold on;
    plot3(lrgpoints(1,:),lrgpoints(2,:),lrgpoints(3,:),'r.','markersize',1);
    plot3(points(1,:),points(2,:),points(3,:),'.','markersize',1);
    plot3(samp(1,:),samp(2,:),samp(3,:),'-*k');
    title(sprintf('Score: %f', err));
    view(20,15);
end


end


function core = rollbackcore(points,lrgpoints,maxvoxdim,rotate)
% something didn't seem to work. Revert back to the first one which
% at least always produces a result..

global nn

[v, ~] = eig(cov(points'));

for k = size(v,2):-1:1
    angle(k) = atan2(norm(cross(v(:,k),[0 0 1])),dot(v(:,k),[0 0 1]));
end
[~, k] = min(angle);

cen = mean(points,2);


% NOTE : theoretically useful to uncomment this if the flame is wide and
% pointing upwards
% I want to put the upwards most vector in the most dominant position (last)
% % % upmost = v(:,k);
% % % v(:,k) = [];
% % % v = [v upmost];

if(~rotate)
    v = eye(3);
end

qrpoints = v'*points;
qrpoints(3,:) = round(qrpoints(3,:)./nn)*nn;
qlrpoints = v'*lrgpoints;
qlrpoints(3,:) = round(qlrpoints(3,:)./nn)*nn;

[a,~] = count_unique(qrpoints(3,:));
core = zeros(3,length(a));
core(3,:) = a';
core(1,:) = arrayfun(@(x)mean(qlrpoints(1,qlrpoints(3,:)==x)),unique(qrpoints(3,:)));
core(2,:) = arrayfun(@(x)mean(qlrpoints(2,qlrpoints(3,:)==x)),unique(qrpoints(3,:)));

core = v*core;

core = fixtopofcore(points,core);
core = fixbottomofcore(points,core);

% core = trimcore(points,core,maxvoxdim);
end


function [bestcore,err] = splinesmooth(core,points,maxbreaks)
global nn
if(isempty(core))
    bestcore = [];
    err = Inf;
    return
end
% figure(5)
% plot3(points(1,:),points(2,:),points(3,:),'.r','markersize',1)
% hold on
% plot3(core(1,:),core(2,:),core(3,:),'*k')

%% find a spline through these points to smooth out the core
x = core(1, :);
y = core(2, :);
z = core(3, :);

t = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)])';
len = max(t);
t = t./len;

bestscore = Inf;
bestcore = [];

qy = points;
qy(3,:) = round(qy(3,:)./nn).*nn;
num = ((max(qy(3,:))-min(qy(3,:)))./nn) + 1;

for breaks = 1:maxbreaks
    pp = splinefit(t,[x ; y ; z],breaks);
    
    % evaluate fit
    samp = ppval(pp,linspace(0,1,num));
    samp(3,:) = round(samp(3,:)./nn).*nn;
    
    dist2core = pdist2(samp',qy','euclidean','smallest',1);
    
    avgscore = mean(dist2core);
    if(avgscore < bestscore)
        bestscore = avgscore;
        bestcore = pp;
    else
        break
    end
    
end

err = bestscore;
end

function [ cut, index ] = pruneplane( points, centre, maxvoxdim )
%PRUNEPLANE Summary of this function goes here
%   Detailed explanation goes here

% % %hhhhack
% % cut = points;
% % index = true(1,size(points,2));

% % return


ndist = 1.5/maxvoxdim; % the threshold for next-to-ness

data = unique(points','rows')';

% choose the point closest to the centre

dists = vnorm(bsxfun(@minus,data,centre));
[d,ix] = min(dists);
if(d > ndist*2)
    error('Prune:NotAdj','The centre point isn''t right next to an actual point');
end

dists = vnorm(bsxfun(@minus,data,data(:,ix)));

neighbours = find(dists < ndist);

if(isempty(neighbours))
    cut = [];
    index = zeros(size(points));
    return;
end

totalix = ix;
while(true)
    totalix = unique([totalix neighbours]);
    ix = neighbours;
    
    rest = data;
    rest(:,totalix) = [];
    restix = 1:size(data,2);
    restix(totalix) = [];
    
    [D,I] = pdist2(rest', data(:,ix)','euclidean','Smallest',8);
    
    neighbours = unique(restix(I(D < ndist)));
    
    if(isempty(neighbours))
        data = data(:,totalix);
        break;
    end
end

if(nargout > 1)
    index = ismember(points',data','rows');
    cut = points(:,index);
else
    cut = points(:,ismember(points',data','rows'));
end

end

function core = fixtopofcore(points,core)
global nn

qpoints = points;
qpoints(3,:) = round(qpoints(3,:)./nn).*nn;
samp = core;
samp(3,:) = round(samp(3,:)./nn).*nn;


for i = size(core,2):-1:1
    d(i) = mean(pdist2(qpoints(:,qpoints(3,:) == samp(3,i))',samp(:,i)','euclidean','Largest',10));
end

topz = max(qpoints(3,:));
if(sum(unique(qpoints(3,:))==topz) < 3)
    topz = topz-nn;
end
% p = mean(qpoints(:,qpoints(3,:) >= topz-2*nn),2);
% dis = vnorm(bsxfun(@minus,samp,p));
dis = topz-samp(3,:);

l = d > dis;
l(1:find(diff(l),1,'last')) = 0;

core(:,l) = [];

end

function core = fixbottomofcore(points,core)
    [dist2core, point2core] = pdist2(core',points','euclidean','smallest',1);
    
    rem = ~ismember(1:length(core),point2core);
    first = find(diff(rem),1,'first');
    if(~isempty(first))
        rem(first+1:end) = 0;
    end
    
    core(:,rem) = [];
end

function core = trimcore(points,core,maxvoxdim)

ndist = 1.5/maxvoxdim; % the threshold for next-to-ness

data = unique(points','rows')';

[D,~] = pdist2(data', core','euclidean','Smallest',1);

core(:,D > ndist) = [];


end