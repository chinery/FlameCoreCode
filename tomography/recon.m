%RECON Reconstruct a 3D voxel density field from input images
% based on Image-Based Tomographic Reconstruction of Flames by I.Ihrke and
% M. Magnor
%
% Here is a basic overview:
% 1. load in the camera calibrations and the images etc
% 2. set up a linear system Sa = p where
%   i. p (nx1) is every pixel from the input images (hence n = num of images * resolution)
%   ii. S (nxm m = total number of voxels) describes the contribution of each
%   voxel to each pixel.
% 3. solve the linear system for a (1xm) which is the density in each voxel
%
% nb: if you are following the paper, then we are using the 'box' basis
% function, so in effect it is one basis function per voxel.
% 
% Computing the matrix S takes a long time and a lot of memory. 
% But you do only have to do it once. 

% The script is set up to be run in parallel.
% Use:
%   matlabpool open x (where x is the number of threads to run (= number of cores))
% before running the file to take advantage of a slight speed increase.
% (assuming you have multiple cores/CPUS)


clear all;
addpath('xml_toolbox','triangle_th','linecurvature_version1b','regu');

%% Set the parameters for the reconstruction.
% best results are supposedly at 64 vox per dim, and I guess the highest
% resolution of image that will fit into memory.
%
% in a later paper the authors update their method to effectively provide much higher
% voxel resolution, which in my own tests do make a big difference. Their update
% used an octree structure which is more efficient in space but requires
% recomputing for each frame. The method implemented here is far less space
% efficient, however has the advantage of only needing to be computed once. I
% have managed to use a resolution of 160 voxels per dimension with 8GB of RAM,
% with the entire original image preserved (ratio = 1).
%
% This script saves the result of S, assuming you are using the GrOVis_Fire
% dataset and calibration, you can make use of one of those. (See below)
imageRatio = 0.4;
voxPerDim = 200;
numOfVoxels = voxPerDim^3;

% where is the data stored
% dataRoot = 'E:/experiment2/data2/';
dataRoot = '/Volumes/CHINERYDATA/experiment2/data2/';
dataFormat = 'png';
numOfCameras = 6;

% Set the boundaries of the voxel space. Should be based on the calibration data. 
global voxelfrom
global voxelto
%voxelfrom = [-0.5,-0.5,-0.5];
% voxelfrom = [0.15,0.0,-0.16];
% voxelfromto = [-1.7 1.3 -0.8 2.2 -1.9 1.1];
%[-1.5 1.5 -1.3 1.7 -1.9 1.1]
% voxelfromto = [-1 1 -0.8 1.2 -1.5 0.5]; % candle
voxelfromto = [-1.5 1.5 -1.3 1.7 -1.5 1.5]; % lighter
%voxelto = [0.5,0.5,0.5];
% voxelto = [0.50,0.47,0.35];
voxelfrom = voxelfromto([1 3 5]);
voxelto = voxelfromto([2 4 6]);

% rotates the cameras
globalRot = [0.908322451910951,0.086804367655184,-0.409164178674602;0.086804367655184,0.917809775663967,0.387414787342509;0.409164178674602,-0.387414787342509,0.826132227574918];

%% read calibration data from xml file
% this is how the data from the paper was stored. Can be replaced with code
% that sets cam(i).K,.R, and .t, corresponding to internal parameters 
% rotation and translation respectively. 
% % % % % % xml = xml_parseany(fscanf(fopen([dataRoot 'tmp/ec.xml']),'%c'));
% % % % % % for i = numOfCameras:-1:1
% % % % % %     tmp = cell2mat(textscan(xml.cam{1,i}.K{1,1}.CONTENT,'%n'));
% % % % % %     cam(i).K = reshape(tmp(3:end),tmp(1),tmp(2))';
% % % % % %     tmp = cell2mat(textscan(xml.cam{1,i}.R{1,1}.CONTENT,'%n'));
% % % % % %     cam(i).R = reshape(tmp(3:end),tmp(1),tmp(2))';
% % % % % %     tmp = cell2mat(textscan(xml.cam{1,i}.t{1,1}.CONTENT,'%n'));
% % % % % %     cam(i).t = reshape(tmp(2:end),tmp(1),1);
% % % % % % end
% % % % % % 
% % % % % % clear xml;
% % % % % % clear tmp;

calibRoot = [dataRoot 'calib/'];
load([calibRoot 'Pmatrices']);
load([calibRoot 'Ce']);
load([calibRoot 'Re']);
for i = numOfCameras:-1:1
    cam(i).P = P((i-1)*3+1:(i-1)*3+3,:);
    cam(i).R = R((i-1)*3+1:(i-1)*3+3,:);
    cam(i).K = cam(i).P(:,1:3) * cam(i).R';
    cam(i).cen = C(:,i);
end

clear P C R

%% added 15/08/2012 to get all fire images...
% this is going to be very platform specific...?
camstr = dataRoot;
camdir = dir([camstr 'cam*']);
for i = length(camdir):-1:1
    imgdirs{i} = dir([camstr camdir(i).name '/*.' dataFormat]);
    ls(i) = length(imgdirs{i});
end
[~, ix] = min(ls);
% % imgdir = imgdirs{ix};


%%
for fireix = 10:length(imgdirs{ix})
%%
fprintf('Reconstructing frame %i',fireix);

%% for loop over colour channels
% alternatively comment out this and the end line to use one channel
% (search for "for loop over colour")
% also set col to something... or uncomment the 'convert to grey' bit just below
% this (and comment out single rgb channel)
for col = 3:-1:1

 clear simages 

for i = length(camdir):-1:1
    colimage = double(imread([camstr camdir(i).name '/' imgdirs{i}(fireix).name ]))./255;
%     % convert to grey
%     images(:,:,i) = rgb2gray(image);
    % single rgb channel
    channel = colimage(:,:,col);
    images(:,:,i) = channel;
end
imgres = size(channel);
clear colimage channel

global S

% % % %% normalise images?
% % % im1 = images(:,:,1);
% % % n = sum(im1(:));
% % % for i = size(images,3):-1:2
% % %     im = images(:,:,i);
% % %     s = sum(im(:));
% % %     images(:,:,i) = (images(:,:,i)./s).*n;
% % % end

%% down-scale images
if imageRatio < 1
    for i = size(images,3):-1:1
        simages(:,:,i) = imresize(images(:,:,i),imageRatio,'lanczos3');
    end
    simages(simages<0) = 0;
    simages(simages>1) = 1;
    
else
    simages = images;
end
clear images;

actualratio = [size(simages,1) size(simages,2)]./imgres;

% %%
% for i = 1:size(simages,3)
%     figure(i);clf;
%     imshow(simages(:,:,i));
% end

%%
% I think the most voxels a single ray can intersect is 4*(voxsize/2)
% but in reality most won't hit the max
maxIntersect = 4*(voxPerDim/2);


%%
% a ray going directly through the diagonal of a single voxel will travel
% this distance
[xmin, xmax, ymin, ymax, zmin, zmax] = getVoxel(1,voxPerDim);
maxVoxSize = norm([xmax,ymax,zmax]-[xmin,ymin,zmin]);



%%
pixels = reshape(simages,[],1);

%%
to = voxelto;
from = voxelfrom;


% colours = hsv(numOfCameras); % this was part of testing the calibration, which I have
% figure; hold on;             % left in because I think it will be useful later

%% load or compute S
% saved .mat files of S are unique per imageRatio and voxPerDim but should
% be cleared out if you are changing the camera calibration or voxel bounds
if exist('S','var') && ~isempty(S)
    % do nothing, this happens in the colour channel loop
elseif exist(sprintf('Sr%3.1fv%2i.mat',imageRatio, voxPerDim),'file')
    load(sprintf('Sr%3.1fv%2i.mat',imageRatio, voxPerDim))
else
    saveram = true;
    if(~saveram)
        for i = size(simages,3):-1:1
            ispars(i).array = zeros(1,size(simages,1)*size(simages,2)*maxIntersect/2);
            jspars(i).array = zeros(1,size(simages,1)*size(simages,2)*maxIntersect/2);
            sspars(i).array = zeros(1,size(simages,1)*size(simages,2)*maxIntersect/2);
            count = zeros(1,size(simages,3));
            %         edges(:,:,i) = edge(simages(:,:,i),'log'); % part of calibration test mentioned earlier
        end
    elseif(saveram && matlabpool('size') > 0)
        error('saveram doesn''t work in parallel');
    else
        try 
            test = spalloc(numel(simages),numOfVoxels,numel(simages)*ceil(maxIntersect/5));
            clear test
        catch Ex
            idSegLast = regexp(Ex.identifier, '(?<=:)\w+$', 'match');
            %clean up for possible error conditions here. Rethrow if unknown error.
            switch idSegLast{1}
                case 'nomem'
                    error('Trying to save RAM, but this size of S just won''t fit');
                otherwise
                    rethrow(Ex)    
            end
        end
        
    end
    
%     edges = reshape(edges,[],1); % part of calibration test
    
% % % %     % supporting data to test coverage
% % % %     vsize = (to-from)./voxPerDim;
% % % %     [x,y,z] = meshgrid(from(1):vsize(1):to(1)-vsize(1),from(2):vsize(2):to(2)-vsize(2),from(3):vsize(3):to(3)-vsize(3));
% % % %     x = x + vsize(1)/2;
% % % %     y = y + vsize(2)/2;
% % % %     z = z + vsize(3)/2;
% % % %     centres = [x(:)'; y(:)'; z(:)'];
% % % %     centres(4,:) = ones(1,length(centres));
% % % % 
% % % %     [y,x] = meshgrid(1:size(simages,1),1:size(simages,2));
% % % %     coverage = [x(:)' ; y(:)'];

    %% ray trace each pixel of each image and calculate distance travelled in each voxel
    % change below line between parfor and for to enable/disable parallel
    for i = 1:size(simages,3) % for each imag
        if(saveram)
            ispars = zeros(1,size(simages,1)*size(simages,2)*maxIntersect/2);
            jspars = zeros(1,size(simages,1)*size(simages,2)*maxIntersect/2);
            sspars = zeros(1,size(simages,1)*size(simages,2)*maxIntersect/2);
            count = 0;
        end
        
%         if(exist(cam(i).P,'var')) % change to is field or sommat
            P = cam(i).P;
            centre = cam(i).cen;
%         else
%             P = cam(i).K * cam(i).R * [ eye(3) -cam(i).t ];
%             centre = cam(i).R * cam(i).t;
%         end
% % % % % % % 
% % % % % % %         %% test which pixels are covered by the voxel grid
% % % % % % %         imgpoints = P*centres;
% % % % % % %         imgpoints(1,:) = imgpoints(1,:)./imgpoints(3,:);
% % % % % % %         imgpoints(2,:) = imgpoints(2,:)./imgpoints(3,:);
% % % % % % %         imgpoints(3,:) = [];
% % % % % % %         imgpoints(:) = imgpoints(:)*imageRatio;
% % % % % % %         
% % % % % % %         imgpoints = floor(imgpoints);
% % % % % % %         uimgpoints = unique(imgpoints','rows')';
% % % % % % %         if(size(imgpoints,2) ~= size(uimgpoints,2))
% % % % % % %             warning('multiple voxels cover one pixel')
% % % % % % %         end
% % % % % % %         clear imgpoints;
% % % % % % %         
% % % % % % %         [~, covered] = ismember(coverage',uimgpoints','rows');
% % % % % % %         covered = reshape(covered,size(simages,1),size(simages,2));
        
        for k = 1:size(simages,2) % for each column of image
            for j = 1:size(simages,1) % for each pixel in column (row)
%                 xy = [(k-1)/imageRatio; (j-1)/imageRatio];
                xy = [(k-1)/actualratio(2); (j-1)/actualratio(1)];
                
                if(saveram)
                    pixelcount = j+(k-1)*size(simages,1);
                else
                    pixelcount = j+(k-1)*size(simages,1)+(i-1)*size(simages,2)*size(simages,1);
                end
                
                % displays status if not run in parallel
                if(mod(pixelcount,5) == 0)
                    clc;
                    fprintf('%f',((j+(k-1)*size(simages,1)+(i-1)*size(simages,2)*size(simages,1))/numel(simages))*100);
                end
                
% % % % % % % %                 if(~covered(k,j))
% % % % % % % %                     continue;
% % % % % % % %                 end
                
                pa = centre;
                apb = P\[xy; 1];
                apb(1:4) = apb(1:4)./apb(4);
                pbv(1:3) = apb(1:3)-centre;
                pbv = pbv/norm(pbv);
                if(norm(centre' + pbv) > norm(centre')) % sometimes the projection moves away from the origin. HACK THAT.
                    pbv = pbv*-1;
                end
                pb = centre' + 200*pbv;
                
                line = (globalRot*[pa'; pb]')';
%                 if(i < 5 && edges(pixelcount) > 0.1) % final part of calibration test
%                     pa = line(1,:);
%                     pb = line(2,:);
%                     plot3([pa(1,1) pb(1,1)], [pa(1,2) pb(1,2)], [pa(1,3) pb(1,3)],'Color',colours(i,:));
%                     plot3(pa(1,1), pa(1,2), pa(1,3),'x','Color',colours(i,:));
%                 end

                %% given a line work out intersection with each of the 
                % planes in each dimension

                % planes = num voxels per dim + 1
                la = line(1,:);
                lb = line(2,:);
                c = 0;
                intersections = zeros(maxIntersect,2);
                interval = (to-from)./voxPerDim;

                for ii = 1:3 % for x, y, z
                    % need three points that are on the plane which we can vary over dim ii
                    base1 = from;
                    if(ii == 1)
                        base2 = from + interval(2)*[0 1 0];
                        base3 = from + interval(3)*[0 0 1];
                    elseif (ii == 2)
                        base2 = from + interval(1)*[1 0 0];
                        base3 = from + interval(3)*[0 0 1];
                    else
                        base2 = from + interval(1)*[1 0 0];
                        base3 = from + interval(2)*[0 1 0];
                    end

                    for jj = 0:voxPerDim % for each plane in this dimension
                        direction = [0 0 0];
                        direction(ii) = 1;
                        % three points on this plane
                        p1 = base1 + jj*interval(ii)*direction;
                        p2 = base2 + jj*interval(ii)*direction;
                        p3 = base3 + jj*interval(ii)*direction;

                        % work out intersections on this plane, add scale factors to list
                        % along with dimension to make working out voxel easier

                        %scale is [t; u; v]
                        scale = [la'-lb', p2'-p1', p3'-p1']\(la'-p1');

                        point = la + scale(1)*(lb-la); %la + t*(lb-la) is the point of intersection
                        if(any(point < from) || any(point > to))
                            continue;
                        end

                        c = c+1;
                        intersections(c,1) = scale(1);
                        intersections(c,2) = ii;
                    end 

                end

                if(any(intersections(:) ~= 0))
                    
                    intersections(intersections(2,:) == 0,:) = [];

                    % sort list of scale factors, go through in pairs finding corresponding
                    % voxel and distance between points and store in S.
                    intersections = sortrows(intersections,1);

                    for ii = 1:(size(intersections,1)-1)
                        p1 = la + intersections(ii,1)*(lb-la);
                        p2 = la + intersections(ii+1,1)*(lb-la);

                        % get unique voxel points for x,y,z
                        if(intersections(ii,2) ~= 1)
                            px = p1(1);
                        else
                            px = p2(1);
                        end
                        if(intersections(ii,2) ~= 2)
                            py = p1(2);
                        else
                            py = p2(2);
                        end
                        if(intersections(ii,2) ~= 3)
                            pz = p1(3);
                        else
                            pz = p2(3);
                        end

                        %indexed by 0
                        xbin = min(floor((px-from(1))/interval(1)),voxPerDim-1);
                        ybin = min(floor((py-from(2))/interval(2)),voxPerDim-1);
                        zbin = min(floor((pz-from(3))/interval(3)),voxPerDim-1);

                        voxNum = voxPerDim^2*zbin + voxPerDim*ybin + xbin + 1;

                        d = norm(p1-p2);
                        if(d ~= 0)
                            if(~saveram)
                                count(i) = count(i) + 1;
                                ispars(i).array(count(i)) = pixelcount;
                                jspars(i).array(count(i)) = voxNum;
                                sspars(i).array(count(i)) = d;
                            else
                                count = count + 1;
                                ispars(count) = pixelcount;
                                jspars(count) = voxNum;
                                sspars(count) = d;
                            end
                        end

                    end
                end
            end
        end
        
        if(saveram)
            if(~exist('C:/tmp','dir'))
                mkdir('C:/tmp');
            end
            save(sprintf('C:/tmp/ispars%i.mat',i),'ispars');
            save(sprintf('C:/tmp/jspars%i.mat',i),'jspars');
            save(sprintf('C:/tmp/sspars%i.mat',i),'sspars');
            clear ispars jspars sspars
        end
    end

    clear count
    
    if(~saveram)
        parfor i = 1:size(simages,3)
            ispars(i).array = ispars(i).array(ispars(i).array ~= 0);
            jspars(i).array = jspars(i).array(jspars(i).array ~= 0);
            sspars(i).array = sspars(i).array(sspars(i).array ~= 0);
        end

        %%
        isparsei = horzcat(ispars(1:size(simages,3)).array);
        clear ispars
        jsparsej = horzcat(jspars(1:size(simages,3)).array);
        clear jspars
        ssparses = horzcat(sspars(1:size(simages,3)).array);
        clear sspars
        S = sparse(isparsei,jsparsej,ssparses,numel(simages),numOfVoxels);
        clear isparsei jsparsej ssparses
    else
        for i = 1:size(simages,3)
            load(sprintf('C:/tmp/ispars%i.mat',i));
            ispars(ispars==0) = [];
            load(sprintf('C:/tmp/jspars%i.mat',i));
            jspars(jspars==0) = [];
            load(sprintf('C:/tmp/sspars%i.mat',i));
            sspars(sspars==0) = [];
            
            Si = sparse(ispars,jspars,sspars,numel(simages)/size(simages,3),numOfVoxels);
            save(sprintf('C:/tmp/Si%i.mat',i),'Si','-v7.3');
            clear Si
        end
        
        clear ispars jspars sspars
        
        S = [];
        for i = 1:6;
            load(sprintf('C:/tmp/Si%i.mat',i));
            S = [S ; Si];
            clear Si
        end
    end
        
    try
        save(sprintf('Sr%3.1fv%2i.mat',imageRatio, voxPerDim),'S','-v7.3')
        if(saveram)
            rmdir('C:/tmp','s');
        end
    catch
    end
end %if exist

%% threshold based on triangle algorithm
% so we can calculate the visual hull a la second part of paper
% camthresh = [0.022, 
for i = size(simages,3):-1:1
    [lehisto ~]=imhist(simages(:,:,i));
    [level]=triangle_th(lehisto,256);
%     threshimages(:,:,i) = imfill(im2bw(simages(:,:,i),max(0,level)),'holes');
    threshimages(:,:,i) = bwareaopen(im2bw(simages(:,:,i),max(0,level)),30);
% % % % % % %     
% % % % % % %     yyy = find(any(threshimages(:,:,i),1));
% % % % % % %     for jj = yyy
% % % % % % %         xxx1 = find(any(threshimages(:,jj,i),2),1,'first');
% % % % % % %         xxx2 = find(any(threshimages(:,jj,i),2),1,'last');
% % % % % % %         threshimages(xxx1:xxx2,jj,i) = 1;
% % % % % % %     end
    
    threshimages(:,:,i) = 1-threshimages(:,:,i);
    
end
thresh = find(reshape(threshimages,[],1));

%%
del = find(any(S(thresh,:),1));

oldS = S;
%% trying to delete the del columns from S without running out of memory
%del = []; %really crap hack to disable visual hull temporarily
% try
    S(:,del) = [];
% catch ME
%     idSegLast = regexp(ME.identifier, '(?<=:)\w+$', 'match');
%     %clean up for possible error conditions here. Rethrow if unknown error.
%     switch idSegLast{1}
%         case 'nomem'
%             S = oldS;
%             for i = size(del,2):-1:1
%                 S(:,del(i)) = [];
%             end
%         otherwise
%             rethrow(ME)    
%     end
% end
    
if(isempty(S))
    cglsix(col) = 0;
    testImage(:,:,col) = zeros(size(simages,1),size(simages,2));
    continue;
end

%% solve the linear system
% several of the inputs and outputs here are unused, see cgls code for more
% info
xitr = 8;%10;
[x,resNE,k,info,allx,lcurve] = cgls('getS',0,pixels,size(S,1),size(S,2),xitr,0.01,1,[]);
% x is the result
% 
lc = lcurve(:,1)./lcurve(:,2);
dlc = lc(2:end)-lc(1:end-1);
ddlc = dlc(2:end)-dlc(1:end-1);
[~,ix] = max(ddlc);

x = allx(:,ix);
% %%
% x = allx(:,1);

% note from 2014: the paper says to use the 'lcurve' to decide which solution to
% pick. I found that ineffective in practice, and never solved why. This part of
% the process was never very intensive, so I think I modified the cgls code to
% return several iterations, then we check below which one most closely matches
% the input visually.

% put the empty voxel columns back into S
S = oldS;
% load(sprintf('Sr%3.1fv%2i.mat',imageRatio, voxPerDim))

% %% check which x gives the best result...
% % % % % % % % % % % % % % % % % for i = xitr:-1:1
% % % % % % % % % % % % % % % % % x = allx(:,i);

% put the empty voxel columns back into a
a = ones(numOfVoxels,1);
a(del) = 0;
a(a == 1) = x;

% 
% % generates images from original positions using a and S
% % newsize = size(simages,1)*size(simages,2);
% % for j = numOfCameras-1:-1:0
% %     %i = 1;
% %     newS = S((j*newsize)+1:(j+1)*newsize,:);
% %     newImage(:,:,j+1) = reshape(newS*a,size(simages,1),[]);
% %     %newImage = reshape(newS*a,size(simages,1),[]);
% % end
% 
% atest = S*a;
% % error(i) = sqrt(mean((simages(:) - atest(:)).^2));
% 
% 
% test = reshape(atest,size(simages,1),size(simages,2),size(simages,3));
% dxtest = convn(test,[-0.5 0 0.5],'same');
% dytest = convn(test,[-0.5 0 0.5]','same');
% mag = (dxtest.^2 + dytest.^2).^(0.5);
% mag = sum(mag(:));
% error(i) = sqrt(mean((simages(:) - atest(:)).^2)) + 0.000002*mag;
% 
% alla(:,i) = a;

% % % % % % % % % % % % % % % % % % % vox = reshape(a,voxPerDim,voxPerDim,voxPerDim);
% % % % % % % % % % % % % % % % % % % colourvox{i}(:,:,:,col) = vox;
% % % % % % % % % % % % % % % % % % % end

%%
% [~, ix] = min(error);
% ix = 4;
% a = alla(:,ix);
fprintf('(actually using %i)',ix);
cglsix(col) = ix;

proj = S*a;
testImage(:,:,col) = reshape(proj(1:size(simages,1)*size(simages,2)),size(simages,1),size(simages,2));

% newsize = size(simages,1)*size(simages,2);
% i = 1; % camera for test
% newS = S((i*newsize)+1:(i+1)*newsize,:);
% testImage(:,:,col) = reshape(newS*a,size(simages,1),[]);

% for xyz = 10:-1:1
%     thisa = alla(:,xyz);
%     testtest{xyz}(:,:,col) = reshape(newS*thisa,size(simages,1),[]);
% end


 
% %% shows the images in the correct order
% for i = 1:8
%     figure; imshow(newImage(:,:,i));
%     %imsave 
% end
 
% %%
% save(['newImage' int2str(col) '.mat'], 'newImage');
% saves the output images. If you are doing all 3 channels, then run
% consolidate_image.m to create the colour images from the saved files

% here is the final result in the correct form!
vox = reshape(a,voxPerDim,voxPerDim,voxPerDim);
colourvox(:,:,:,col) = vox;

%% show a 3D voxel representation using vol3d
% figure;vol3d('cdata',vox)

end %of for loop over colour channels

% for xyz = 10:-1:1
%     figure(xyz);clf;
%     imshow(testtest{xyz});
% %     imwrite(testtest{xyz},sprintf('testiteration%i.png',xyz));
% end
% return
% figure;vol3d('cdata',colourvox./10,'alpha',sum(colourvox./10,4)./3);
% return;

if(~exist('testImage','var'))
    continue;
end
imwrite(testImage,sprintf([dataRoot 'volumes/test/frame%.5ivox%i-%i-%i.png'],fireix,cglsix(1),cglsix(2),cglsix(3)),'png')
% imwrite(testImage,sprintf('frame%.3ivox%i.png',fireix,ix),'png')
% the below added 'colourvox' to a file already containing vox
% load(sprintf('../AllVolumes/frame%.3ivox.mat',fireix));
% save(sprintf('../AllVolumes/frame%.3ivox.mat',fireix),'vox', 'colourvox');

% assuming this script is always run in colour mode and no grey version is ever
% needed, this is more useful
% nb: make sure the directory exists, and/or change it to one you like.
% save(sprintf('../AllVolumes/frame%.3ivox.mat',fireix), 'colourvox');
save(sprintf([dataRoot 'volumes/frame%.5ivox.mat'],fireix), 'colourvox','-v7.3');


end % of loop over all images