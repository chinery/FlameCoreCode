data_root = '../GenerateFlames2014/results/animations/cbatch19/animation1/smooth/img/';
% data_root = '../MyData2014/lighter1/smoothflames/img/proc/';
% data_root = '../ImageMerge/vid6/source4/';
images = dir([data_root '*.png']);

output = [data_root 'proc/'];
% output = 'C:/Users/Andrew/Dropbox/Work/PhD/ImageMerge/vid6/source4/';
if(~exist(output,'dir'))
    mkdir(output);
end

rs = [901 993];
% rs = [412 549];

% % crop = [412 549 273 160]; %height, width, y, x
% crop = [700 933 1 1];
crop = [507 676 145 167];
cropfirst = false;

% crop = [789 1052 51 88];
% cropfirst = true;

rs2 = [480 640];

for i = 1:length(images)
    im = imread([data_root images(i).name]);
    
    if(cropfirst)
        im = im(crop(3):crop(3)+crop(1)-1,crop(4):crop(4)+crop(2)-1,:);
    end
    
    im = imresize(im,rs);
    
    if(~cropfirst)
        im = im(crop(3):crop(3)+crop(1)-1,crop(4):crop(4)+crop(2)-1,:);
    end
    
    if(~isempty(rs2))
        im = imresize(im,rs2);
    end
    
    imwrite(im,[output images(i).name]);
end
