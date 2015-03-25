%% smooth
clear all;
% data_root = '../MyData2014/candle1-2/';
% data_root = '../GenerateFlames2014/results/animations/cbatch21/animation1/';
% data_root = '../../MyData2014/lighter1/';
% data_root = '/Volumes/CHINERYDATA/experiment2/data1/';
data_root = '../generation/results/batch22/animation1/';
% data_root = '../data/fitting/twod/';
% data_root = 'F:/experiment2/data1/';

% resmode: work out of a folder (e.g. generated) rather than on source data
% (lots of folders)
resmode = true; 

if(~resmode)
    flames = dir([data_root 'flames/*.mat']);
else
    flames = dir([data_root '*.mat']);
end
limit = length(flames);
%  limit = 30;
weight = ones(1,limit);
crud = false(1,length(flames));
crud([3819, 3827:3831, 3836:3845, 3919:3923, 4226:4230, 4279:4377, 4401:4424, 4551:4575, 4593:4663, 5017:5029, 5247:5264, 5342:5355]) = true;
for frame = limit:-1:1

    if(~resmode)
        load([data_root 'flames/' flames(frame).name]);
    else
        load([data_root flames(frame).name]);
    end
    if(isfield(flame,'normm'))
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
    else
        for i = 3:-1:1
            allnormm{i}(:,frame) = zeros(50,1);
            allnorms{i}(:,frame) = zeros(50,1);
            allunifo{i}(:,frame) = zeros(50,1);
            allscale{i}(:,frame) = zeros(50,1);
            allmaxim{i}(:,frame) = zeros(50,1);
        end
        weight(frame) = 0;
    end
    if(crud(frame))
        weight(frame) = 0;
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

save all.mat
figure
%%
for frame = 1:limit
    normm = [normmfit{1}(:,frame)'; normmfit{2}(:,frame)'; normmfit{3}(:,frame)'];
    norms = [normsfit{1}(:,frame)'; normsfit{2}(:,frame)'; normsfit{3}(:,frame)'];
    unifo = [unifofit{1}(:,frame)'; unifofit{2}(:,frame)'; unifofit{3}(:,frame)'];
    scale = [scalefit{1}(:,frame)'; scalefit{2}(:,frame)'; scalefit{3}(:,frame)'];
    maxim = [maximfit{1}(:,frame)'; maximfit{2}(:,frame)'; maximfit{3}(:,frame)'];
    if(~resmode)
        load([data_root 'flames/' flames(frame).name]);
    else
        load([data_root flames(frame).name]);
    end
    flame.normm{1} = normm;
    flame.norms{1} = norms;
    flame.unifo{1} = unifo;
    flame.scale{1} = scale;
    flame.maxim{1} = maxim;
    if(crud(frame) || ~isfield(flame,'cores'))
        flame.discard(1) = true;
    end
    if(~resmode)
        save(sprintf([data_root 'smoothflames/frame%05i.mat'],frame),'flame')
    else
        save(sprintf([data_root 'smooth/frame%05i.mat'],frame),'flame')
    end
    if(crud(frame) || ~isfield(flame,'cores'))
        continue
    end
    displaysingleflame(flame);view(80,0);camproj('perspective');
    view(80,0);
    if(~resmode)
        print('-dpng',sprintf([data_root 'smoothflames/img/frame%05i'],frame),'-opengl');
    else
        print('-dpng',sprintf([data_root 'smooth/img/frame%05i'],frame),'-opengl');
    end
end
% 
% generateflames
% %%
% clear flames
% for frame = 10:-1:1
%     load(sprintf([data_root 'nmpismoothflames/frame%05i.mat'],frame),'flame')
%     flames{frame} = flame;
% end
% figure(1);clf;
% displaysingleflame(flames,[data_root 'nmsmoothflames/tests/frame%05itesst.png']);
