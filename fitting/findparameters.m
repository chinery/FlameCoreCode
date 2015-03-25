% once the cores have been found, this code will find the other parameters

run('../setup.m');
addpath([pwd '/helpers/']);

% data_root = '../MyData2014/candle1-2/';
% data_root = '../MyData2014/lighter1/';
data_root = 'F:/experiment2/data1/';
% data_root = '/Volumes/CHINERYDATA/experiment2/data1/';
data_root = '/Users/andy/copy/work/PhD/MyData2014/lighter1/';

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

% this lets you renumber the flames, no matter what you set here, the first one will be number 1 when saved
startix = 1; 

global iii
iii = 1;
for frame = startix:length(cores)
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
        defival = 20;   %  <---- the default number of points along the core to use
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
        ival = defival;  % <-- was once based on length (above), now just set

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
% %     lastflame = flame;
% % %     
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

% now you can run smoothflames to smooth the parameters if necessary
