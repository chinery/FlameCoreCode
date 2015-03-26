clearvars -except data_root

% batch convert from vox to points
% data_root = '../MyData2014/candle1-2/';
% data_root = 'E:/experiment2/data1/';
vol_root = [data_root 'volumes/'];
volumes = dir([vol_root '*.mat']);

for frame = 10:length(volumes)
    fprintf('Converting frame %i\n',frame);
    load(sprintf([vol_root 'frame%05ivox.mat'],frame));
    
    colourvox = colourvox./10;
    
    vox = sum(colourvox,4)./3;
    points = vox2points(vox,1,100,0);
    
% %     figure;
% %     plot3(points(1,:),points(2,:),points(3,:),'.','markersize',1);
% %     view(-41,-5);
% %     axis([0.300000000000000,0.600000000000000,0.400000000000000,0.700000000000000,0.150000000000000,0.650000000000000;])
% %     return
    
    spoints = vox2points(vox,1,10);
    for i = 3:-1:1
        colourpoints{i} = vox2points(colourvox(:,:,:,i),1,100);
    end
    
    save(sprintf([data_root 'points/frame%05ipoints.mat'],frame),'points','spoints','colourpoints');
end