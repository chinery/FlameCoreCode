frame = 207;
newnumber = 2;

if(exist('numberoflicks','var'))
    clear numberoflicks;
end
load(sprintf('../AllVolumes/frame%.3ivox.mat',frame));

numberoflicks = newnumber;

save(sprintf('../AllVolumes/frame%.3ivox.mat',frame),'colourvox','vox','numberoflicks');