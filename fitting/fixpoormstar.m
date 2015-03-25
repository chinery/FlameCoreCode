% fix poor mstar

data_root = 'F:/experiment2/data1/';
% data_root = '/Volumes/CHINERYDATA/experiment2/data1/';
flames = dir([data_root 'smoothtrimflames/*.mat']);
limit = length(flames);
% limit = 4474;


%%
for frame = limit:-1:1
%     load([data_root 'flames/' flames(frame).name]);
    load([data_root 'smoothtrimflames/' flames(frame).name]);
%     load([data_root flames(frame).name]);
    if(isfield(flame,'normm') && isfield(flame,'cores'))
        mstar = (flame.maxim{1}-flame.normm{1})./flame.norms{1};
        if(mstar(1,end)>mstar(1,end-1))
            mstar(:,end) = mstar(:,end-1);
            flame.maxim{1} = mstar.*flame.norms{1} + flame.normm{1};
       
            fprintf('Frame %i FIXED\n',frame);
        else
            fprintf('Frame %i no fix\n',frame);
        end
%         displaysingleflame(flame);view(80,0);camproj('perspective');
%         print('-dpng',sprintf([data_root 'smoothtrimfixflames/img/frame%05i'],frame),'-opengl');
    else
        fprintf('Frame %i no flame\n',frame);
    end
    save([data_root 'smoothtrimfixflames/' flames(frame).name],'flame')
end    