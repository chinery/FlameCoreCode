clear all

% flames_root = '../MyData2014/candle1-2/nmcgcsmoothflames/';
% flames_root = '../MyData2014/lighter1/smoothflames/';
flames_root = 'F:/experiment2/data1/smoothtrimflames/';
% flames_root = '/CHINERYDATA/experiment2/data1/smoothflames';
files = dir([flames_root '*.mat']);

% control = load('pinkflame','pinkflame');
% control = control.pinkflame;
control = [];

points = [];
scalepoints = [];
rotpoints = []; 
transpoints = [];

normmpoints = [];
normspoints = [];
unifopoints = [];
densipoints = [];
maxespoints = [];
match = [];

mapback = [];

defival = 10;
rsize = @(x)(interp1(x',linspace(1,size(x,2),defival))');

% refsamp = bsxfun(@minus,refsamp,mean(refsamp,2));
% refsamp = bsxfun(@minus,refsamp,refsamp(:,1)); 
conseq = false;
c = 1;
for i = 1:length(files)
    load([flames_root files(i).name]);
       
    if(~isfield(flame,'cores'))
        conseq = false;
        continue;
    end
    for j = 1:length(flame.cores)
%             if(c > length(filt))
%                 return;
%             end
%             if(filt(c))
%                 c = c + 1;
%                 continue;
%             end
%             c = c + 1;
            if(flame.discard(j))
                conseq = false;
                continue;
            end
            if(conseq)
                match(1,end+1) = c-1;
                match(2,end) = c;
            end
            conseq = true;
            c = c + 1;
            
            samp = flame.cores{j};
            if(size(samp,2) ~= defival)
                samp = interp1(samp',linspace(1,size(samp,2),defival),'spline')';
            end

            newsamp = samp;
            
            mv = [eye(3) , -samp(1:3,1)];
            mv(4,:) = [0 0 0 1];
            R = [1 0 0 0;...
                 0 1 0 0;...
                 0 0 1 0;...
                 0 0 0 1];
            mvb = [eye(3) , samp(1:3,1)];
            mvb(4,:) = [0 0 0 1];

            T = mvb*R*mv;
            
            newsamp(4,:) = 1;
            newsamp = T*newsamp;
            newsamp = bsxfun(@rdivide, newsamp, newsamp(4,:));
            newsamp = newsamp(1:3,:);
            
            points = [points, reshape(newsamp,[],1)];
%             scalepoints = [scalepoints, scale'];
%             rotpoints = [rotpoints, rotation'];
%             transpoints = [transpoints, mx];

            normm = reshape(rsize(flame.normm{j}(1,:)),[],1);
            norms = reshape(rsize(flame.norms{j}(1,:)),[],1);
            unifo = reshape(rsize(flame.unifo{j}),[],1);
            densi = reshape(rsize(flame.scale{j}),[],1);
            maxes = reshape(rsize(flame.maxim{j}(1,:)),[],1);
            
%             normm = reshape(rsize(flame.normm{j}(1,:)),[],1)+normmmap;
%             norms = reshape(rsize(flame.norms{j}(1,:)),[],1)+normsmap;
%             unifo = reshape(rsize(flame.unifo{j}),[],1)+unifomap;
%             densi = reshape(rsize(flame.scale{j}),[],1)+densimap;
%             maxes = reshape(rsize(flame.maxim{j}(1,:)),[],1)+maxesmap;
            
%             if(i==9)
%                 sum(reshape(control.normm{1}(1,:),[],1)-normm)
%                 sum(reshape(control.norms{1}(1,:),[],1)-norms)
%                 sum(reshape(control.maxim{1}(1,:),[],1)-maxes)
%                 sum(reshape(control.unifo{1},[],1)-unifo)
%                 sum(reshape(control.scale{1},[],1)-densi)
%                 pause
%             end
            
            normmpoints = [normmpoints, normm];
            normspoints = [normspoints, norms];
            unifopoints = [unifopoints, unifo];
            densipoints = [densipoints, densi];
            maxespoints = [maxespoints, maxes];
            
% %         end
    end
    
end

if(~isempty(control))
    ctrlpoint = reshape(rsize(control.scale{1}),[],1);
%     [~,trainingstartix] = pdist2(densipoints',ctrlpoint','euclidean','Smallest',1);
    [~,trainingstartix] = pdist2(densipoints',zeros(size(ctrlpoint))','euclidean','Smallest',1);
% trainingstartix = 9;
    densimap = ctrlpoint - densipoints(:,trainingstartix);
    
%     ctrlpoint = reshape(control.cores{1},[],1);
%     [~,trainingstartix] = pdist2(points',ctrlpoint','euclidean','Smallest',1);
%     startix = trainingstartix;
    
%     ctrlpoint = rsize(control.normm{1}(1,:));
%     [~,startix] = pdist2(normmpoints',ctrlpoint','euclidean','Smallest',1);
%     normmmap = ctrlpoint - normmpoints(:,startix);
%     
%     ctrlpoint = rsize(control.norms{1}(1,:));
%     [~,startix] = pdist2(normspoints',ctrlpoint','euclidean','Smallest',1);
%     normsmap = ctrlpoint - normspoints(:,startix);
%     
% %     maxespoints = (maxespoints - normmpoints)./normspoints;
%     
%     ctrlpoint = (rsize(control.maxim{1}(1,:))-rsize(control.normm{1}(1,:)))./rsize(control.norms{1}(1,:));
%     [~,startix] = pdist2(maxespoints',ctrlpoint','euclidean','Smallest',1);
%     maxesmap = ctrlpoint - maxespoints(:,startix);
    
    ctrlpoint = reshape(rsize(control.unifo{1}),[],1);
%     [~,startix] = pdist2(unifopoints',ctrlpoint','euclidean','Smallest',1);
    [~,startix] = pdist2(unifopoints',zeros(size(ctrlpoint))','euclidean','Smallest',1);
    unifomap = ctrlpoint - unifopoints(:,startix);
    

    
%     normmpoints = max(bsxfun(@plus,normmpoints,normmmap),0);
%     normspoints = max(bsxfun(@plus,normspoints,normsmap),0);
%     maxespoints = max(bsxfun(@plus,maxespoints,maxesmap),0);
    unifopoints = max(bsxfun(@plus,unifopoints,unifomap),0);
    densipoints = max(bsxfun(@plus,densipoints,densimap),0);
    
%     maxespoints = maxespoints.*normspoints + normmpoints;
    
    block = normspoints <= 0 | maxespoints <= 0 | reshape(all(reshape(densipoints,3,[]) == 0,1) | any(reshape(densipoints,3,[]) < 0,1),defival,[]);
    densipoints(reshape(repmat(reshape(block,1,[]),3,1),defival*3,[])) = 0;
    maxespoints(block) = 0;
    
end

% match = [1:size(points,2)-1; 2:size(points,2)];
% match = [1:size(points,2)-4; 5:size(points,2)]; %25fps

% make matching 25fps
tempmatch = match;
for i = 1:3
    [lia,ix] = ismember(tempmatch(2,:),match(1,:));
    tempmatch = [tempmatch(1,lia); match(2,ix(ix>0))];
end
match = tempmatch;

% save ncolourtraining25fps points normmpoints normspoints unifopoints densipoints maxespoints match
save latestcandle10p25fps points normmpoints normspoints unifopoints densipoints maxespoints match
% save pinklightertraining points normmpoints normspoints unifopoints densipoints maxespoints match trainingstartix
% save colourtrainingtall points normmpoints normspoints unifopoints densipoints maxespoints match
