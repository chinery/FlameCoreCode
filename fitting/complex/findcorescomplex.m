clear all

% This is to find cores in complex flames (more than one lick per flame)
% once findparts has already been run
% please see notes in findcores.m in parent folder

useoldmethod = true;
rotate = true;

% data_root = '../MyData2014/candle1-2/';
% data_root = '../MyData2014/lighter1/';
% data_root = '/Volumes/CHINERYDATA/experiment2/data1/';
% data_root = 'G:/experiment2/data1/';
data_root = '/Users/andy/copy/work/PhD/MyData2014/alcohol1/';
pointsdir = dir([data_root 'points/*.mat']);
%%
for frame = 31:length(pointsdir)
    frame
    clearvars -except frame data_root pointsdir
    
    load(sprintf([data_root 'volumes/' pointsdir(frame).name],frame));
    
    if(numberoflicks == 0)
        continue;
    end
    flame.voxsize = size(vox);
    clear vox colourvox
    
    load([data_root 'points/' pointsdir(frame).name]);
    
    load([data_root 'parts/' sprintf('frame%03i.mat',frame)]);
%     ix = ones(1,length(points));
    
    success = 0;
    discarded = [];
    flame.discard = zeros(1,max(ix));
    for i = 1:max(ix)-1
        part = points(:,ix==i);
        
        part = cleanuppoints( part, flame.voxsize(1) );
        
        if(all(part(3,:) > 0.5))
            discard = 1;
            flat = 0;
        else
            discard = 0;
            nn = 1/flame.voxsize(1);
%             [x,y,z] = meshgrid(min(part(1,:)):nn:max(part(1,:)),min(part(2,:)):nn:max(part(2,:)),(min(part(3,:))-5*nn):nn:(min(part(3,:)+10*nn)));
%             lrgpart = [part [x(:)';y(:)';z(:)']];%part;%points(:,ix==i | ix==max(ix));
            
%             [core,flat] = fitcore(part,[],flame.voxsize(1));
            
            lrgpart = points(:,ix==i | ix==max(ix));
            [core,flat] = fitcore(part,lrgpart,flame.voxsize(1),true);
        end
        
        flat = 0;
        if(flat)
            figure(1);clf;hold on;
            plot3(part(1,:),part(2,:),part(3,:),'.r');
            samp = ppval(core,linspace(0,1,100));
            plot3(samp(1,:),samp(2,:),samp(3,:),'*k');
            
            def = { num2str(1) };
            options.WindowStyle = 'normal';
            x = inputdlg('Discard? \n 1: Yes 0: Keep', 'Bad Sector?', 1, def,options);
            if(isempty(x))
                return
            end
            discard = str2double(x{:});
            close(1)
        end
        
        flame.discard(i) = discard;
        
        if(~discard && ~isempty(core))
            success = success+1;
            flame.cores{success} = core;
            flame.flat(success) = false; % kinda redundant... 
        else
            discarded = [discarded i];
        end
            
    end
    
    flame.discard(max(ix)) = true;
% %     
    figure(1);clf;
    hold on;
    colours = hsv(max(ix));
    j = 0;
    for i = 1:max(ix)-1
        if(ismember(i,discarded))
            continue;
        end
        j = j + 1;
        part = points(:,ix==i);
        plot3(part(1,:),part(2,:),part(3,:),'.','color',colours(i,:),'markersize',1);
        samp = ppval(flame.cores{j},linspace(0,1,100));
%         samp = flame.cores{j};
        plot3(samp(1,:),samp(2,:),samp(3,:),'*k');
    end
    view(20,15);
%     axis([0.300000000000000,0.600000000000000,0.400000000000000,0.700000000000000,0.150000000000000,0.650000000000000]);
    axis([0 1 0 1 0 1]);

    % f = splinesmooth, added bottom of core fix, improved top
    print('-dpng',sprintf([data_root 'othercores/img/frame%05if'],frame));
    save(sprintf([data_root 'othercores/frame%05i.mat'],frame),'flame')
    
end

% % %%
% % cores = dir([data_root 'cores/*.mat']);
% % for frame = 1:length(cores)
% %     if(frame == 1)
% %         load(sprintf([data_root 'cores/frame%05i.mat'],frame));
% %         lastcore = ppval(flame.cores{1},linspace(0,1,100));
% %         load(sprintf([data_root 'cores/frame%05i.mat'],frame+1));
% %         nextcore = ppval(flame.cores{1},linspace(0,1,100));
% %     else
% %         core = nextcore;
% %         load(sprintf([data_root 'cores/frame%05i.mat'],frame+1));
% %         nextcore = ppval(flame.cores{1},linspace(0,1,100));
% %         
% %         core = (1/3).*lastcore + (1/3).*core + (1/3).*nextcore;
% %         
% %         load(sprintf([data_root 'points/frame%05ipoints.mat'],frame));
% %         
% %         figure(1);clf;hold on;
% %         part = points;
% %         plot3(part(1,:),part(2,:),part(3,:),'.','color',[1 0 0],'markersize',1);
% %         plot3(core(1,:),core(2,:),core(3,:),'*k');
% %         view(20,15);
% %         axis([0.300000000000000,0.600000000000000,0.400000000000000,0.700000000000000,0.150000000000000,0.650000000000000]);
% %         print('-dpng',sprintf([data_root 'cores/img/frame%05it'],frame));
% %         lastcore = core;
% %     end
% % end
return

