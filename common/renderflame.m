function [ vox ] = renderflame( flame, dim, ptsoncore, fix )
%RENDERFLAME this is like newrenderskeleton but for the NEW model
%(normm,norms,unifo,scale)


if nargin < 2
    if ~isfield(flame,'reconvox')
        error('Need to provide the voxel resultion, e.g. [160 160 160]');
    end
    dim = size(flame.reconvox);
end

if length(dim) == 1
    dim = [dim dim dim];
end

if nargin < 3 || isempty(ptsoncore)
    ptsoncore = 50;
end

if nargin < 4
    fix = true;
end

vox = zeros([dim 3]);

[x,y,z] = meshgrid( 0:dim(1), 0:dim(2), 0:dim(3));
voxpoints = bsxfun(@rdivide,[x(:)' ; y(:)' ; z(:)'],dim');

for i = 1:length(flame.cores)
    if(isfield('flame','flat') && flame.flat(i))
        continue;
    end
    
    
    normm = flame.normm{i};
    norms = flame.norms{i};
    scale = flame.scale{i};
    unifo = flame.unifo{i};
    maxim = flame.maxim{i};
    
    maxes = max(flame.maxim{i},[],1);
    
    
    if(fix)
        % max in terms of std away from mean
        mstar = (maxim-normm)./norms;
        if(mstar(1,end)>mstar(1,end-1))
            mstar(:,end) = mstar(:,end-1);
            maxim = mstar.*norms + normm;
        end
    end
    
    if(length(scale) ~= ptsoncore)
        scale = interp1(scale',linspace(1,length(scale),ptsoncore))';
        normm = interp1(normm',linspace(1,length(normm),ptsoncore))';
        norms = interp1(norms',linspace(1,length(norms),ptsoncore))';
        unifo = interp1(unifo',linspace(1,length(unifo),ptsoncore))';
        maxim = interp1(maxim',linspace(1,length(maxim),ptsoncore))';

        maxes = interp1(maxes',linspace(1,length(maxes),ptsoncore))';
    end
    
        

    
    if(isstruct(flame.cores{i}))
        sampcore = ppval(flame.cores{i},linspace(0,1,ptsoncore));
    else
        sampcore = flame.cores{i};
        if(size(sampcore,2) ~= ptsoncore)
            sampcore = interp1(sampcore',linspace(0,size(sampcore,2),ptsoncore),'spline','extrap')';
        end
        x = sampcore(1, :); 
        y = sampcore(2, :); 
        z = sampcore(3, :);

        t = cumsum([0;sqrt(diff(x(:)).^2 + diff(y(:)).^2 + diff(z(:)).^2)])';
        len = max(t);
        t = t./len;
        splinecore = splinefit(t,[x ; y ; z],4);
    end
    
    
    sampcore(:,end) = sampcore(:,end)-0.9*(sampcore(:,end)-sampcore(:,end-1));
    
    if(~isfield(flame,'flat'))
        nsc = sampcore;
    elseif(flame.flat(i))
        nsc = sampcore;
    else
        nsc = sampcore(:,1)-0.1*(sampcore(:,2)-sampcore(:,1));
        nsc(:,2:ptsoncore+1) = sampcore;
%         nsc(:,end+1) = sampcore(:,end)+0.1*(sampcore(:,end)-sampcore(:,end-1));
    end
    
    
%     figure; hold on; view(20,15);
    %% go along core and find bbox
    if(isstruct(flame.cores{i}))
        cprime = ppdiff(flame.cores{i});
        corevecs = ppval(cprime,linspace(0,1,ptsoncore));
        cprimeprime = ppdiff(cprime);
        corenorms = ppval(cprimeprime,linspace(0,1,ptsoncore));
    else
        cprime = ppdiff(splinecore);
        corevecs = ppval(cprime,linspace(0,1,ptsoncore));
        cprimeprime = ppdiff(cprime);
        corenorms = ppval(cprimeprime,linspace(0,1,ptsoncore));
    end
    
    normcorevecs = bsxfun(@rdivide,corevecs,vnorm(corevecs));
    vecdists = [];
    vecdists(2:ptsoncore) = vnorm(sampcore(:,2:end) - sampcore(:,1:end-1));
    vecdists(1) = vecdists(2);
    
    normcorenorms = bsxfun(@rdivide,corenorms,vnorm(corenorms));
    
    bbox = [Inf 0 ; Inf 0 ; Inf 0];
    
    
    for ci = 1:ptsoncore
        tan = normcorevecs(:,ci);
        nor = normcorenorms(:,ci);
        bin = cross(nor,tan);
        if(norm(bin) ~= 0 && ~any(isnan(bin)))
            bin = bin./norm(bin);
        else
            continue;
        end
        
        fframe = [bin nor tan];
        
        [x, y] = pol2cart(linspace(0,2*pi,100),maxes(ci));
        rvecs = fframe'\[x; y; zeros(1,100)];
        addpoints = bsxfun(@plus,sampcore(:,ci), bsxfun(@times,rand(1,100)*vecdists(ci),corevecs(:,ci))) + rvecs;
        if(ci == ptsoncore) % allow for the 'ploom' at the end
            ffframe = [tan nor bin];
            rrvecs = ffframe'\[x; y; zeros(1,100)];
            addpoints = [addpoints bsxfun(@plus,sampcore(:,ci), bsxfun(@times,rand(1,100)*vecdists(ci),corevecs(:,ci))) + rrvecs];
        end
        tbb = [min(addpoints,[],2), max(addpoints,[],2)];
        bbox = [min(bbox(:,1),tbb(:,1)) max(bbox(:,2),tbb(:,2))];
    end
    
    %% remove points outside of bbox
    bblog = (voxpoints(1,:) < bbox(1,1) ...
           | voxpoints(2,:) < bbox(2,1) ...
           | voxpoints(3,:) < bbox(3,1) ...
           | voxpoints(1,:) > bbox(1,2) ...
           | voxpoints(2,:) > bbox(2,2) ...
           | voxpoints(3,:) > bbox(3,2));
       
    cutvp = voxpoints(:,~bblog);
            
    
    %% calculate closest points
    [dis2cor, p2c] = pdist2(nsc',cutvp','euclidean','smallest',1);
    
    point2core = ones(size(voxpoints,2),1);
    dist2core = Inf(size(voxpoints,2),1);
    point2core(~bblog) = p2c;
    dist2core(~bblog) = dis2cor;
    
    if(isfield(flame,'flat') && ~flame.flat(i))
        delete = point2core == 1;
        point2core = max(point2core - 1,1);
        point2core = min(point2core,ptsoncore);
    end
    
    step = @(m,u,x)(heaviside(m-x).*u);
    
    modelfun = @(mu,sig,u,m,s,x)(s.*(...
        (step(m,1-u,x)).*(...
            exp(-((x-mu).^2)./(2*sig.^2)))...
        + step(m,u,x))); 
    
    
    for col = 3:-1:1
        thisvox = modelfun(normm(col,point2core)',norms(col,point2core)',unifo(col,point2core)',maxim(col,point2core)',scale(col,point2core)',dist2core);

        if(isfield(flame,'flat') && ~flame.flat(i))
            thisvox(delete) = 0;
        end
        thisvox(isnan(thisvox)) = 0;

        thisvox = reshape(thisvox,dim(1)+1,dim(2)+1,dim(3)+1);
        [x, y, z] = meshgrid(linspace(1,dim(1)+1,dim(1)),linspace(1,dim(2)+1,dim(2)),linspace(1,dim(3)+1,dim(3)));
        thisvox = interp3(thisvox,x,y,z);

        corevox(:,:,:,col) = thisvox;
    end
    
    vox = vox + corevox;
end


end

