% once you have run numberlicks, this code will try to find the splits
% this code hasn't been used in a long time!! so might need some structural
% editing
% but you should get the gist from what's here.
% also I never had great source data for complex flames, so I don't guarantee
% it'll work WELL.

data_root = '/Users/andy/copy/work/PhD/MyData2014/alcohol1/';

pointdir = dir([data_root 'points/*.mat']);

for frame = 1:length(pointdir)
    frame
    clearvars -except frame framestart data_root pointdir
%     if(exist(sprintf('../AllVolumes/points/parts/frame%.3i.mat',frame),'file'))
%         continue;
%     end
    
    if(exist('numberoflicks','var'))
        clear numberoflicks;
    end
    
    load(sprintf([data_root 'volumes/' pointdir(frame).name],frame));
    
    if(~exist('numberoflicks','var'))
        error('need to run numberlicks for this frame first');
    end
    
    if(numberoflicks == 0)
        continue;
    end
    
        load(sprintf([data_root 'points/' pointdir(frame).name],frame));
    
   
    % k means still sucks.
    % ix = kmeans(points',num+1);

    % %%
    % mix = GMM_fit( points, num+1, 100, 1e-3 );
    % ix = GMM_classify( points, mix );

    % %%
    % ix = emgm(points,num+1);

    %%
    scale = 15;

    p2 = points;
    cutf = 0.2;
    p2(:,p2(3,:) < cutf) = [];
    
    siz = size(vox,1);
    p2 = cleanuppoints(p2,siz);
    
    oldp2 = p2;

    p2(3,:) = p2(3,:)./scale;
        
    figure(1);clf;
    vol3d('cdata',min(colourvox*6,1),'alpha',min(sum(colourvox,4)*2,1));
    title(sprintf('Number of licks: %i',numberoflicks));
    view(20,15);
    
    smallsamp = spoints;
    smallsamp(:,smallsamp(3,:) < cutf) = [];
    guessnum = guesslicks( smallsamp, siz, false );
    if(guessnum == numberoflicks)
%         ix = emgm(p2,numberoflicks);
        [~,ix] = chincluster(oldp2,size(vox,1));
    elseif(guessnum == 1)
%          ix = emgm(p2,numberoflicks);
        [~,ix] = definedmeanshiftcluster(p2,numberoflicks);
    else
        [~,chinix] = chincluster(oldp2,size(vox,1));
        while(true)
            if(max(chinix) >= numberoflicks)
                if(numberoflicks < max(chinix))
                    numberoflicks = max(chinix);
                end
                ix = chinix;
                break;
            end
            h = figure;
            hold on;
            colours = hsv(max(chinix));
            dia = 'Which should be split?';
            for i = 1:max(chinix)
                part = oldp2(:,chinix==i);
                plot3(part(1,:),part(2,:),part(3,:),'.','color',colours(i,:));

                width(i) = norm(max(part(1:2,:),[],2)-min(part(1:2,:),[],2));

                dia = [dia sprintf('\n%i:-  r: %i g: %i b: %i',i,colours(i,1),colours(i,2),colours(i,3))];
            end
            view(20,15);

            [~,mix] = max(width);
            def = { num2str(mix) };
            options.WindowStyle = 'normal';
            options.Interpreter='tex';
            x = inputdlg(dia, 'Split', 1, def,options);
            if(isempty(x))
                return
            end
            choose = str2double(x{:});


            part = oldp2(:,chinix==choose);

            if(guessnum == numberoflicks - 1)
                close(h);

                [~,ix] = definedmeanshiftcluster(part,2);

                ix(ix==2) = max(chinix) + 1;
                ix(ix==1) = choose;
                chinix(chinix == choose) = ix;

                ix = chinix;
                
                if(numberoflicks < max(chinix))
                    numberoflicks = max(chinix);
                end
                
                break;
            else
                def = { num2str(2) };
                options.WindowStyle = 'normal';
                options.Interpreter='tex';
                split = Inf;
                while(split > numberoflicks - max(chinix) + 1)
                    x = inputdlg('How many parts should it be split into?', 'Split', 1, def,options);
                    if(isempty(x))
                        return
                    end
                    split = str2double(x{:});
                end
                
                [~,ix] = definedmeanshiftcluster(part,split);
                
                % chinix has 1s and 2s
                % number of licks is 5
                % chinix 1 is broken into 3
                % so ix goes from 1 to 3
                % want to set ix to 1,3,4
                
                newlabels = [choose, max(chinix)+1:max(chinix)+split-1];
                for j = split:-1:1
                    ix(ix==j) = newlabels(j);
                end
                chinix(chinix == choose) = ix;
            end
        end
        
        
    end
    
    [~,loc] = ismember(points',oldp2','rows');
    loc(loc~=0) = ix(loc(loc~=0));
    loc(loc==0) = numberoflicks+1;
    init = loc';

    backup = points(3,:);
    points(3,:) = points(3,:)./scale;
    ix = emgm(points,init);
% ix = init;
    points(3,:) = backup;
    


    figure(2);clf;
    hold on;
    colours = hsv(numberoflicks+1);
%     colours = [ 1 0 1; 0 1 1];
    for i = 1:numberoflicks
        part = points(:,ix==i);
        plot3(part(1,:),part(2,:),part(3,:),'.','color',colours(i,:),'markersize',1);
    end
    view(20,15);

    pause;
    save([data_root sprintf('/points/parts/frame%.3i.mat',frame)],'ix');
end


