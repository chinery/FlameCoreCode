% this is based on newfitparameters which is the new model but now being more
% smart with colour stuff

run('../setup.m');
addpath([pwd '/helpers/']);

% data_root = '../MyData2014/candle1-2/';
% data_root = '../MyData2014/lighter1/';
data_root = 'F:/experiment2/data1/';
% data_root = '/Volumes/CHINERYDATA/experiment2/data1/';

cores = dir([data_root 'smoothcores/*.mat']);
pointdir = dir([data_root 'points/*.mat']);

% finds all the parameters for each core.
% allgams = [];
% allweib = [];
% allgaus = [];
% load(sprintf([data_root 'nmpiflames/frame%05i.mat'],333),'flame')
% lastflame = flame;
if(~exist('lastflame','var'))
    lastflame = [];
end
if(~exist('knn','var'))
    global knn
    load('fitclassifier');
end

startix = 1; % this is to renumber the flames
global iii
iii = 1;
for frame = startix+5333:length(cores)
%     try
    frame
    clearvars -except frame data_root lastflame knn cores pointdir startix iii%allgams allweib allgaus
    
%     load(sprintf('../AllVolumes/frame%.3ivox.mat',frame),'numberoflicks');
%     
%     if(numberoflicks == 0)
%         continue;
%     end
    numberoflicks = 1;
    
    load([data_root 'points/' pointdir(frame).name]);
    
%     load(sprintf('../AllVolumes/points/parts/frame%.3i.mat',frame));
%     clusterix = ix;
    clusterix = ones(1,length(points));
    
    load([data_root 'smoothcores/' cores(frame).name]);
    
    if(~isfield(flame,'cores') || isempty(flame.cores{1}))
        continue;
    end
    
    clustrange = 1:max(clusterix);
%     clustrange(logical(flame.discard)) = [];
    
    maxvoxdim = max(flame.voxsize);
    for i = 1:length(clustrange)
        j = clustrange(i);
        greypart = points(:,clusterix == j);
        % clean up dodgy data
%         if(length(clustrange) > 1)
%             greypart = cleanuppoints( greypart, flame.voxsize(1) );
%         end

        % prepare core
        defival = 50;
        if(isstruct(flame.cores{i}))
            sampcore = ppval(flame.cores{i},linspace(0,1,defival));
        else
            sampcore = flame.cores{i};
        end
% % % 
% % %             len = sum(vnorm(sampcore(:,2:end) - sampcore(:,1:end-1)));
% % %             ival = min(ceil(len/(1/maxvoxdim)),200);
% % %             ival = max(10,ival);
            
% % %             ival = round(ival/2);
        ival = defival;

        if(isstruct(flame.cores{i}))
            sampcore = ppval(flame.cores{i},linspace(0,1,ival));
        else
            sampcore = interp1(sampcore',linspace(0,size(sampcore,2),ival),'spline','extrap')';
        end
        
        part = greypart;
        
        [dist2core, point2core] = pdist2(sampcore',part','euclidean','smallest',1);
        vec2core = vnorm(part - sampcore(:,point2core));

        %%
        point2core(dist2core > (mean(dist2core) + 2*std(dist2core))) = [];
        vec2core(:,dist2core > (mean(dist2core) + 2*std(dist2core))) = [];
        cutpart = part;
        cutpart(:,dist2core > (mean(dist2core) + 2*std(dist2core))) = [];
        dist2core(dist2core > (mean(dist2core) + 2*std(dist2core))) = [];

        %%
        [c,ia,ic] = unique(cutpart','rows');
        % c' = part(:,ia)
        % part = c'(:,ic)
        uniquep2c = point2core(ia);
        [a,b] = count_unique(uniquep2c);
        voxperpoint = zeros(length(sampcore),1);
        voxperpoint(a) = b;

        %% fit parameters
        [a,b] = count_unique(ic);
        % the number of points at cutpart(:,i) is b(i)
        sf = 100; %scale factor used when converting to points
        dens = b./sf; 

        uniqdists = dist2core(ia);
        uniquep2c = point2core(ia);

        range = find(voxperpoint > 20);
        
                    
%             figure(1);clf;
%             hist(dist2core,20);
%             export_fig(sprintf('./results/hists2/%03i.png',iii),'-transparent')
%             iii = iii + 1;
            
%             arrayfun(@(x)histandsave(dist2core(point2core==x)),range);
%             continue
        
        paramfit = arrayfun(@(x)estimateparameters(uniqdists(uniquep2c==x)',dens(uniquep2c==x)),range);
        
        a = range;
        normm = zeros(1,length(sampcore));
        norms = zeros(1,length(sampcore));
        unifo = zeros(1,length(sampcore));
        scale = zeros(1,length(sampcore));
        maxim = zeros(1,length(sampcore));

        normm(a) = arrayfun(@(x)x.mu,paramfit);
        norms(a) = arrayfun(@(x)x.sig,paramfit);
        unifo(a) = arrayfun(@(x)x.u,paramfit);
        scale(a) = arrayfun(@(x)x.s,paramfit);
        maxim(a) = arrayfun(@(x)x.m,paramfit);
        
        flame.normm{i}(1:3,:) = repmat(normm,3,1);
        flame.norms{i}(1:3,:) = repmat(norms,3,1);
        flame.maxim{i}(1:3,:) = repmat(maxim,3,1);
        
        for col = 3:-1:1

            mem = ismember(colourpoints{col}',greypart','rows');
            part = colourpoints{col}(:,mem);
            
            % work out the closest point on the core for each point in the part
            
            [dist2core, point2core] = pdist2(sampcore',part','euclidean','smallest',1);
            vec2core = vnorm(part - sampcore(:,point2core));

            %%
            point2core(dist2core > (mean(dist2core) + 2*std(dist2core))) = [];
            vec2core(:,dist2core > (mean(dist2core) + 2*std(dist2core))) = [];
            cutpart = part;
            cutpart(:,dist2core > (mean(dist2core) + 2*std(dist2core))) = [];
            dist2core(dist2core > (mean(dist2core) + 2*std(dist2core))) = [];
            
            %%
            [c,ia,ic] = unique(cutpart','rows');
            % c' = part(:,ia)
            % part = c'(:,ic)
            uniquep2c = point2core(ia);
            [a,b] = count_unique(uniquep2c);
            voxperpoint = zeros(length(sampcore),1);
            voxperpoint(a) = b;

            %% fit parameters
            [a,b] = count_unique(ic);
            % the number of points at cutpart(:,i) is b(i)
            sf = 100; %scale factor used when converting to points
            dens = b./sf; 
            
            uniqdists = dist2core(ia);
            
            range = find(voxperpoint > 20);
            
            scalefit = arrayfun(@(x)fitscaleparameter(uniqdists(uniquep2c==x)',dens(uniquep2c==x),...
                normm(x),norms(x),unifo(x),scale(x),maxim(x)),range);
            
            a = range;
            cunifo = zeros(1,length(sampcore));
            cscale = zeros(1,length(sampcore));
            cunifo(a) = arrayfun(@(x)x.u,scalefit);
            cscale(a) = arrayfun(@(x)x.s,scalefit);
            
            flame.unifo{i}(col,:) = cunifo;
            flame.scale{i}(col,:) = cscale;
            
            % hoping to remove this later
            maxes = zeros(length(sampcore),1);
            maxes(a) = arrayfun(@(x)max(dist2core(point2core==x)),range);
            
            flame.maxes{i}(col,:) = maxes;
        end
        
%         end % of col loop

    end   


    save(sprintf([data_root 'flames/frame%05i.mat'],frame-startix+1),'flame')
%     lastflame = flame;
% %     
    figure(1);clf;
    displaysingleflame(flame);
    camproj('perspective')
    view(80,0);
    saveas(gcf,sprintf([data_root 'flames/img/frame%05i.png'],frame-startix+1));
%     fprintf('save disabled!')
    
%     catch errobj
%         fprintf('Error: %s', getReport(errobj));
%         emailChineryError(errobj);
%     end

end 

% generateflames

% 
% [ v2  ] = newrenderskeleton( flame, flame.voxsize(1), 300, true );
% 
% figure;
% vol3d('cdata',v2,'alpha',sum(v2,4)./3);view(20,15)

return
%% smooth
clear all;
data_root = '../MyData2014/candle1-2/';
flames = dir([data_root 'nmcgcflames/*.mat']);
limit = length(flames);
weight = ones(1,limit);
for frame = limit:-1:1
    load(sprintf([data_root 'nmcgcflames/frame%05i.mat'],frame));
    for i = 3:-1:1
        allnormm{i}(:,frame) = flame.normm{1}(i,:)';
        allnorms{i}(:,frame) = flame.norms{1}(i,:)';
        allunifo{i}(:,frame) = flame.unifo{1}(i,:)';
        allscale{i}(:,frame) = flame.scale{1}(i,:)';
        allmaxim{i}(:,frame) = flame.maxim{1}(i,:)';
    end
    if(mean(flame.scale{1}(1,flame.scale{1}(1,:) ~= 0 & flame.scale{1}(2,:) ~= 0)./flame.scale{1}(2,flame.scale{1}(1,:) ~= 0 & flame.scale{1}(2,:) ~= 0)) < 1.32)
        weight(frame) = 0.2;
    end
end    
    
t1 = 1:size(allnormm{1},1);
% t1 = t1./200;
t1 = repmat(t1',1,limit);
% t1 = t1.*10;

t2 = 1:limit;
% t2 = t2./length(flames);
t2 = repmat(t2,size(allnormm{1},1),1);

normmfit = allnormm;
normsfit = allnorms;
unifofit = allunifo;
scalefit = allscale;
maximfit = allmaxim;

window = 10;
overlap = 9;
for i = 3:-1:1
%     normmfit{i} = multisplinesmooth(t1,allnormm{i});
    normmfit{i} = surfacesmooth(t1,t2,normmfit{i},5,4);
%     normsfit{i} = multisplinesmooth(t1,allnorms{i});
    normsfit{i} = surfacesmooth(t1,t2,normsfit{i},5,4);
%     unifofit{i} = multisplinesmooth(t1,allunifo{i});
    unifofit{i} = surfacesmooth(t1,t2,unifofit{i},5,4);
%     scalefit{i} = multisplinesmooth(t1,allscale{i});
    scalefit{i} = surfacesmooth(t1,t2,scalefit{i},window,overlap,weight);
%     maximfit{i} = multisplinesmooth(t1,allmaxim{i});
    maximfit{i} = surfacesmooth(t1,t2,maximfit{i},5,4);
end

for frame = 1:limit
    normm = [normmfit{1}(:,frame)'; normmfit{2}(:,frame)'; normmfit{3}(:,frame)'];
    norms = [normsfit{1}(:,frame)'; normsfit{2}(:,frame)'; normsfit{3}(:,frame)'];
    unifo = [unifofit{1}(:,frame)'; unifofit{2}(:,frame)'; unifofit{3}(:,frame)'];
    scale = [scalefit{1}(:,frame)'; scalefit{2}(:,frame)'; scalefit{3}(:,frame)'];
    maxim = [maximfit{1}(:,frame)'; maximfit{2}(:,frame)'; maximfit{3}(:,frame)'];
    load(sprintf([data_root 'nmcgcflames/frame%05i.mat'],frame));
    flame.normm{1} = normm;
    flame.norms{1} = norms;
    flame.unifo{1} = unifo;
    flame.scale{1} = scale;
    flame.maxim{1} = maxim;
    save(sprintf([data_root 'nmcgcsmoothflames/frame%05i.mat'],frame),'flame')
    displaysingleflame(flame);view(80,0);camproj('perspective');
    print('-dpng',sprintf([data_root 'nmcgcsmoothflames/img/frame%05inms'],frame),'-opengl');
end
% 
% generateflames
return
%%
clear flames
for frame = 10:-1:1
    load(sprintf([data_root 'nmpismoothflames/frame%05i.mat'],frame),'flame')
    flames{frame} = flame;
end
figure(1);clf;
displaysingleflame(flames,[data_root 'nmsmoothflames/tests/frame%05itesst.png']);
