clear all

run('../setup.m');

load('../data/generation/latestcandle10p.mat');
% load('strongerfilter.mat');
% load('thefilter.mat');
% load('../data/generation/ncolourtraining25fps.mat');
% load('../data/generation/pinklightertraining.mat');


% a logical vector will filter out points
% e.g. using the first 2000 of 'latestcandle10p'
% has less erractic motion
thefilter = false(1,size(points,2));
thefilter(2001:end) = true;



% we are going to roughly follow a 'path' through each parameter space
% with some controlled randomness.
% with all parameters set to 0 (or whatever), this script should start
% at a random flame from the training, and then just follow the original video

% % pick a start point from training
% % look at n nearest training points within distance d
% % pick one at random weighted by closeness
% % (do not average, [1 0] and [0 1] may be valid but [0.5 0.5] not)
% % pick a number of points in a gaussian area around the matching 'next' point
% % % sigma parameter, possibly derrived from distance between data
% % pick one randomly, weighted by density estimate and gaussian liklihood

parameters.animation.totalitr = 1;
parameters.animation.anilength = 2000;
parameters.animation.start = 1; % 0 is random

parameters.animation.dirstring = './results/batch22/animation%i/';
parameters.nearest.num = 30;
parameters.nearest.distance = 1.5;
% parameters.nearest.exactdistcorrection = 0.1;
parameters.nearest.vmfkappa = 10;
parameters.gaussian.num = 20;
parameters.gaussian.scale = 0.3; % 0 to 1 probably

% if lock is enabled, then the next point from training will be found for the
% first parameter, and then the same point will be used in the other parameters
parameters.lock = false;
% advanced kernel density estimation requires the code found at:
% http://www.vicos.si/File:Maggot_v3.4.zip
% which needs to be compiled for your system
parameters.useadvancedkde = false;

parameters.disablesave = false;
if(parameters.disablesave)
    warning('save disabled!');
end

modifier.shape = 1.5;
modifier.normms = 0.4;
modifier.norms = 0.4;
modifier.normm = 0.4;
modifier.maxes = 0.4;
modifier.densi = 0.4;
modifier.unifo = 0.4;
modifier.original = parameters.gaussian.scale;

% this highlight feature never worked perfectly
highlightix = 0; % 0 to disable
% highlightix = 3580;

% % % set up

% filter the data (should be done before saving really)
points(:,thefilter) = [];
normmpoints(:,thefilter) = [];
normspoints(:,thefilter) = [];
maxespoints(:,thefilter) = [];
densipoints(:,thefilter) = [];
unifopoints(:,thefilter) = [];
delmap = 1:length(thefilter);
delmap(thefilter) = [];
del = find(thefilter);
match(:,ismember(match(1,:),del)) = [];
match(:,ismember(match(2,:),del)) = [];
[~,match(1,:)] = ismember(match(1,:),delmap);
[~,match(2,:)] = ismember(match(2,:),delmap);
[~,highlightix] = ismember(highlightix,delmap);
[~,parameters.animation.start] = ismember(parameters.animation.start,delmap);


% shape
eigmodels.shape = EigModel_make(points,'keepf',0.970);
modelpoints.shape = EigModel_project(points,eigmodels.shape);
if(~parameters.useadvancedkde)
    Hshape = densitykernel(modelpoints.shape);
    densityfunction.shape = @(test)densityestimate(test,modelpoints.shape,Hshape);
else
    pdfshape = advanceddensitykernel(modelpoints.shape);
    densityfunction.shape = @(test)advanceddensityestimate(test,pdfshape);
end
if(highlightix > 0)
    shapehighlight.point = modelpoints.shape(:,highlightix);
    shapehighlight.maxdist = pdist2(modelpoints.shape',shapehighlight.point','euclidean','Largest',1);
    nextpointweight.shape = @(test)simpleweight(test,shapehighlight);
end
fprintf('...\n');


% density based parameters,
% normal
eigmodels.normm = EigModel_make( normmpoints, 'keepf', 0.97 );
modelpoints.normm = EigModel_project(normmpoints,eigmodels.normm);
% Hnormm = densitykernel(modelpoints.normm);
% densityfunction.normm = @(test)densityestimate(test,modelpoints.normm,Hnormm);


eigmodels.norms = EigModel_make( normspoints, 'keepf', 0.97 );
modelpoints.norms = EigModel_project(normspoints,eigmodels.norms);
% Hnorms = densitykernel(modelpoints.norms);
% densityfunction.norms = @(test)densityestimate(test,modelpoints.norms,Hnorms);

% maximum
normspoints(normspoints == 0) = 1e10;
maxespoints(isnan(maxespoints)) = 0;
maxespoints = (maxespoints - normmpoints)./normspoints;
eigmodels.maxes = EigModel_make( maxespoints, 'keepf', 0.97 );
modelpoints.maxes = EigModel_project(maxespoints,eigmodels.maxes);
% Hmaxes = densitykernel(modelpoints.maxes);
% densityfunction.maxes = @(test)densityestimate(test,modelpoints.maxes,Hmaxes);

% GROUP the above three
combo = [modelpoints.normm ; modelpoints.norms ; modelpoints.maxes];
whiten = std(combo,0,2);
combo = bsxfun(@rdivide,combo,whiten);

eigmodels.normms = EigModel_make( combo, 'keepf', 1 );
modelpoints.normms = EigModel_project(combo,eigmodels.normms);
if(~parameters.useadvancedkde)
    Hnormms = densitykernel(modelpoints.normms);
    densityfunction.normms = @(test)densityestimate(test,modelpoints.normms,Hnormms);
else
    pdfshape = advanceddensitykernel(modelpoints.normms);
    densityfunction.normms = @(test)advanceddensityestimate(test,pdfshape);
end
fprintf('...\n');

% density (scale)
eigmodels.densi = EigModel_make( densipoints, 'keepf', 0.970 );
modelpoints.densi = EigModel_project(densipoints,eigmodels.densi);
if(~parameters.useadvancedkde)
    Hdensi = densitykernel(modelpoints.densi);
    densityfunction.densi = @(test)densityestimate(test,modelpoints.densi,Hdensi);
else
    pdfshape = advanceddensitykernel(modelpoints.densi);
    densityfunction.densi = @(test)advanceddensityestimate(test,pdfshape);
end
fprintf('...\n');

% uniform ratio
eigmodels.unifo = EigModel_make( unifopoints, 'keepf', 0.970 );
modelpoints.unifo = EigModel_project(unifopoints,eigmodels.unifo);
if(~parameters.useadvancedkde)
    Hunifo = densitykernel(modelpoints.unifo);
    densityfunction.unifo = @(test)densityestimate(test,modelpoints.unifo,Hunifo);
else
    pdfshape = advanceddensitykernel(modelpoints.unifo);
    densityfunction.unifo = @(test)advanceddensityestimate(test,pdfshape);
end
fprintf('...\n');
% HACKY HACK
params = {'shape', 'normms', 'densi', 'unifo'};
% params = {'shape', 'normm', 'norms', 'maxes', 'densi', 'unifo'};

% load strongfiltdensityfn

for i = 1:length(params)
%     eval(sprintf('sigmas.%s = diag(mean((diff(modelpoints.%s,1,2)).^2,2));',params{i},params{i}));
    eval(sprintf('difpoints = modelpoints.%s(:,match(2,:)) - modelpoints.%s(:,match(1,:));',params{i},params{i}));
    eval(sprintf('temp = diag(mean((difpoints).^2,2));'));
    eval(sprintf('sigmas.%s = cov(modelpoints.%s'');',params{i},params{i}));
    eval(sprintf('sigmas.%s = (sigmas.%s./sigmas.%s(1,1)).*temp(1,1);',params{i},params{i},params{i}));
%     eval(sprintf('dists.%s = mean(vnorm(difpoints))+std(vnorm(diff(modelpoints.shape,1,2)));',params{i},params{i}));
    eval(sprintf('dists.%s = mean(vnorm(difpoints));',params{i},params{i}));
end

% % path
% if(highlightix > 0)
%     score = pathtohighlightscore(modelpoints.shape, match, highlightix, dists.shape*2, parameters );
%     if(any(score == -1))
%         badix = find(score == -1);
%         match(:,ismember(match(2,:),badix)) = [];
%     end
%     
%     meanscore = mean(score(score~=-1));
%     worst = max(score)+round(0.05*max(score));
%     score(score == -1) = worst;
%     shapehighlight.scaledscore = score./meanscore;
%     nextpointweight.shape = @(test)pathweight(modelpoints.shape,test,shapehighlight.scaledscore);
% %     shapehighlight.point = modelpoints.shape(:,highlightix);
% %     shapehighlight.maxdist = pdist2(modelpoints.shape',shapehighlight.point','euclidean','Largest',1);
% %     nextpointweight.shape = @(test)simpleweight(test,shapehighlight);
% end

% make matching 25fps
tempmatch = match;
for i = 1:3
    [lia,ix] = ismember(tempmatch(2,:),match(1,:));
    tempmatch = [tempmatch(1,lia); match(2,ix(ix>0))];
end
match = tempmatch;

if(highlightix > 0)
    for i = 2:length(params)
        eval(sprintf('nextpointweight.%s = [];',params{i}));
    end
else
    for i = 1:length(params)
        eval(sprintf('nextpointweight.%s = [];',params{i}));
    end
end

%%
fprintf('starting\n')
% %% make the animation
for iteration = 1:parameters.animation.totalitr
    if(~exist(sprintf(parameters.animation.dirstring,iteration),'dir'))
        mkdir(sprintf(parameters.animation.dirstring,iteration));
    end
    
    if(parameters.animation.start == 0)
        startix = randi(size(modelpoints.shape,2));
    else
        startix = parameters.animation.start;
    end
    
    for i = 1:length(params)
        eval(sprintf('%s = modelpoints.%s(:,startix);',params{i},params{i}));
        eval(sprintf('direction.%s = zeros(size(%s));',params{i},params{i}));
    end
    
    for frame = 1:parameters.animation.anilength
        %% save a frame
        unprojshape = EigModel_unproject(shape,eigmodels.shape);
        newcore = reshape(unprojshape,3,[]);
        
        upnormms = EigModel_unproject(normms,eigmodels.normms)'.*whiten'; 
        upnormm = EigModel_unproject(upnormms(1:size(modelpoints.normm,1))',eigmodels.normm)';
        upnorms = EigModel_unproject(upnormms(size(modelpoints.normm,1)+1:size(modelpoints.normm,1)+size(modelpoints.norms,1))',eigmodels.norms)';
        upmax = (EigModel_unproject(upnormms(size(modelpoints.normm,1)+size(modelpoints.norms,1)+1:size(modelpoints.normm,1)+size(modelpoints.norms,1)+size(modelpoints.maxes,1))',eigmodels.maxes)'.*upnorms)+upnormm;
        
%         upnormm = EigModel_unproject(normm,eigmodels.normm)';
%         upnorms = EigModel_unproject(norms,eigmodels.norms)';
%         upmax = EigModel_unproject(maxes,eigmodels.maxes)';
%         upmax = (upmax.*upnorms)+upnormm;

        upnormm = repmat(upnormm,3,1);
        upnorms = repmat(upnorms,3,1);
        upmax = repmat(upmax,3,1);
        
        updens = reshape(EigModel_unproject(densi,eigmodels.densi),3,[]);
        upunifo = reshape(EigModel_unproject(unifo,eigmodels.unifo),3,[]);
        
        updens = max(0,updens);
        upmax = max(0,upmax);
        
        figure(12);clf;
        flame = displaysinglecore(newcore,upnormm,upnorms,updens,upmax,upunifo);
        
        if(~parameters.disablesave)
            save(sprintf([parameters.animation.dirstring '/%04iparams.mat'],iteration,frame),'flame');
            print('-dpng',sprintf([parameters.animation.dirstring '%04i'],iteration,frame)); 
        else
            fprintf('.');
            if(randi(30,1) == 1)
                fprintf('\n');
            end
            drawnow
        end
        if(frame == parameters.animation.anilength)
            break;
        end
        
        %% get next frame
        for i = 1:length(params)
            eval(sprintf('parameters.gaussian.scale = modifier.original * modifier.%s;',params{i}));
            if(~parameters.lock)
                eval(sprintf('[%s, direction.%s] = getnextpoint(%s,direction.%s,modelpoints.%s,match,sigmas.%s,dists.%s,densityfunction.%s,nextpointweight.%s,parameters);',params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i}))
            elseif(parameters.lock && i == 1)
                eval(sprintf('[%s, direction.%s,lockix] = getnextpoint(%s,direction.%s,modelpoints.%s,match,sigmas.%s,dists.%s,densityfunction.%s,nextpointweight.%s,parameters);',params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i}))
            else
                eval(sprintf('[%s, direction.%s] = getnextpoint(%s,direction.%s,modelpoints.%s,match,sigmas.%s,dists.%s,densityfunction.%s,nextpointweight.%s,parameters,lockix);',params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i}))
            end
        end
    end
end