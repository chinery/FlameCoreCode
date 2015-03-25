data_root = 'F:/experiment2/data1/';
flames = dir([data_root 'flames/*.mat']);
limit = length(flames);

weight = ones(1,limit);

%%
rsize = @(x)(interp1(x',linspace(1,size(x,2),defival))');
for frame = limit:-1:1
    load([data_root 'flames/' flames(frame).name]);
%     load([data_root flames(frame).name]);
    if(isfield(flame,'normm') && isfield(flame,'cores'))
        prob = all(flame.scale{1} < 0.0001,1) | all(flame.norms{1} < 0.0001,1);
        trims = find(~prob,1,'first');
        trime = find(~prob,1,'last');
        
        if(length(flame.cores{1}) ~= size(flame.scale{1},2))
            flame.cores{1} = interp1(flame.cores{1}',linspace(1,length(flame.cores{1}),size(flame.scale{1},2)))';
        end
        
        if(trims ~= 1 || trime ~= size(flame.scale{1},2))
            flame.cores{1} = interp1(flame.cores{1}',linspace(trims,trime,size(flame.cores{1},2)))';
            flame.normm{1} = interp1(flame.normm{1}',linspace(trims,trime,size(flame.normm{1},2)))';
            flame.norms{1} = interp1(flame.norms{1}',linspace(trims,trime,size(flame.norms{1},2)))';
            flame.scale{1} = interp1(flame.scale{1}',linspace(trims,trime,size(flame.scale{1},2)))';
            flame.unifo{1} = interp1(flame.unifo{1}',linspace(trims,trime,size(flame.unifo{1},2)))';
            
            if(trims ~= 1 && trime ~= size(flame.scale{1},2))
                fprintf('Frame %i trim both!\n',frame);
            elseif(trims ~= 1)
                fprintf('Frame %i trim bottom!\n',frame);
            else
                fprintf('Frame %i trim top!\n',frame);
            end
        else
            fprintf('Frame %i no trim\n',frame);
        end
        displaysingleflame(flame);view(80,0);camproj('perspective');
        print('-dpng',sprintf([data_root 'trimflames/img/frame%05i'],frame),'-opengl');
    else
        fprintf('Frame %i no flame\n',frame);
    end
    save([data_root 'trimflames/' flames(frame).name],'flame')
end    