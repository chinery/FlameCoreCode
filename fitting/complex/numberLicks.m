% run this first for a complex (multiple lick) flame
% the code is really old, from when I only had very bad data and was trying
% anything to get better results, including manually labelling the number of
% licks in each frame.
% you probably won't want to do this(!!) especially for 100fps footage
% you might be able to streamline the process

%%
for frame = 30:300;
% frame = 115;

if(exist('numberoflicks','var'))
    clear numberoflicks;
end
load(sprintf('../AllVolumes/frame%.3ivox.mat',frame));

if(~exist('numberoflicks','var'))
    load(sprintf('../AllVolumes/points/frame%.3ivox.mat',frame));
    guess = guesslicks(spoints, size(vox,1));
    
    figure(1);clf;
    vol3d('cdata',min(colourvox*6,1),'alpha',min(sum(colourvox,4)*2,1));
    title(sprintf('Frame: %i',frame));
    view(20,15);
    
    def = { num2str(guess) };
    options.WindowStyle = 'normal';
    x = inputdlg('Number of licks:', 'Licks', 1, def,options);
    if(isempty(x))
        return
    end
    numberoflicks = str2double(x{:});
    
    save(sprintf('../AllVolumes/frame%.3ivox.mat',frame),'colourvox','vox','numberoflicks');
end

end