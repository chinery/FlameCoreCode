% data_root = '../generation/results/batch7/animation1/smooth/';
data_root = '../../MyData2014/rrr/';
% data_root = '../data/fitting/twod/smoothflames/';
out = [data_root 'img/'];

flames = dir([data_root '*.mat']);

angle = 80;
changed = 0;
for i = 1:length(flames)
    load([data_root flames(i).name],'flame')
    
    figure(12);clf;
    if(isfield(flame,'cores') && ~ isempty(flame.cores) &&  ~exist(sprintf([out 'frame%05i.png'],str2num(flames(i).name(6:10))),'file'))
    
        flame.scale{1} = flame.scale{1}*1.2;
        displaysingleflame(flame);
        axis equal
%         if(i == start)
%             a = axis;
%         end
%         axis(a)
        view(80,0);
    %     view(angle-i,0)
    end
    if( ~exist(sprintf([out 'frame%05i.png'],str2num(flames(i).name(6:10))),'file'))
       print('-dpng',sprintf([out 'frame%05i'],str2num([flames(i).name(6:10)]))); 
       changed = changed + 1;
    end
    
end