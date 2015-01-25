clear all

run('../setup.m');

load('../data/generation/ncolourtraining25fps.mat');
% load('../data/generation/pinklightertraining.mat');

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
parameters.animation.start = 220; % 0 is random
parameters.animation.dirstring = './results/batch1/animation%i/';
parameters.nearest.num = 10;
parameters.nearest.distance = 1;
parameters.nearest.exactdistcorrection = 0.1;
parameters.gaussian.num = 20;
parameters.gaussian.scale = 0.3; % 0 to 1 probably

modifier.shape = 1.5;
modifier.normms = 0.4;
modifier.densi = 0.4;
modifier.unifo = 0.4;
modifier.original = parameters.gaussian.scale;

highlightix = 0; % 0 to disable
% highlightix = 220;

% % % set up
% shape
eigmodels.shape = EigModel_make(points,'keepf',0.90);
modelpoints.shape = EigModel_project(points,eigmodels.shape);
Hshape = densitykernel(modelpoints.shape);
densityfunction.shape = @(test)densityestimate(test,modelpoints.shape,Hshape);
if(highlightix > 0)
    shapehighlight.point = modelpoints.shape(:,highlightix);
    shapehighlight.maxdist = pdist2(modelpoints.shape',shapehighlight.point','euclidean','Largest',1);
    nextpointweight.shape = @(test)simpleweight(test,shapehighlight);
end

% density based parameters,
% normal
eigmodels.normm = EigModel_make( normmpoints, 'keepf', 0.9 );
modelpoints.normm = EigModel_project(normmpoints,eigmodels.normm);
Hnormm = densitykernel(modelpoints.normm);
densityfunction.normm = @(test)densityestimate(test,modelpoints.normm,Hnormm);

eigmodels.norms = EigModel_make( normspoints, 'keepf', 0.9 );
modelpoints.norms = EigModel_project(normspoints,eigmodels.norms);
Hnorms = densitykernel(modelpoints.norms);
densityfunction.norms = @(test)densityestimate(test,modelpoints.norms,Hnorms);

% maximum
normspoints(normspoints == 0) = eps;
maxespoints = (maxespoints - normmpoints)./normspoints;

eigmodels.maxes = EigModel_make( maxespoints, 'keepf', 0.9 );
modelpoints.maxes = EigModel_project(maxespoints,eigmodels.maxes);
Hmaxes = densitykernel(modelpoints.maxes);
densityfunction.maxes = @(test)densityestimate(test,modelpoints.maxes,Hmaxes);

% GROUP the above three
combo = [modelpoints.normm ; modelpoints.norms ; modelpoints.maxes];
whiten = std(combo,0,2);
combo = bsxfun(@rdivide,combo,whiten);

eigmodels.normms = EigModel_make( combo, 'keepf', 1 );
modelpoints.normms = EigModel_project(combo,eigmodels.normms);
Hnormms = densitykernel(modelpoints.normms);
densityfunction.normms = @(test)densityestimate(test,modelpoints.normms,Hnormms);

% density (scale)
eigmodels.densi = EigModel_make( densipoints, 'keepf', 0.90 );
modelpoints.densi = EigModel_project(densipoints,eigmodels.densi);
Hdensi = densitykernel(modelpoints.densi);
densityfunction.densi = @(test)densityestimate(test,modelpoints.densi,Hdensi);

% uniform ratio
eigmodels.unifo = EigModel_make( unifopoints, 'keepf', 0.90 );
modelpoints.unifo = EigModel_project(unifopoints,eigmodels.unifo);
Hunifo = densitykernel(modelpoints.unifo);
densityfunction.unifo = @(test)densityestimate(test,modelpoints.unifo,Hunifo);

% HACKY HACK
params = {'shape', 'normms', 'densi', 'unifo'};

for i = 1:length(params)
    eval(sprintf('sigmas.%s = diag(mean((diff(modelpoints.%s,1,2)).^2,2));',params{i},params{i}));
    eval(sprintf('dists.%s = mean(vnorm(diff(modelpoints.%s,1,2)));',params{i},params{i}));
end

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
        
        upnormm = repmat(upnormm,3,1);
        upnorms = repmat(upnorms,3,1);
        upmax = repmat(upmax,3,1);
        
        updens = reshape(EigModel_unproject(densi,eigmodels.densi),3,[]);
        upunifo = reshape(EigModel_unproject(unifo,eigmodels.unifo),3,[]);
        
        updens = max(0,updens);
        upmax = max(0,upmax);
        
        figure(12);clf;
        flame = displaysinglecore(newcore,upnormm,upnorms,updens,upmax,upunifo);
        
        save(sprintf([parameters.animation.dirstring '/%04iparams.mat'],iteration,frame),'flame');
        print('-dpng',sprintf([parameters.animation.dirstring '%04i'],iteration,frame)); 
        if(frame == parameters.animation.anilength)
            break;
        end
        
        %% get next frame
        for i = 1:length(params)
            eval(sprintf('parameters.gaussian.scale = modifier.original * modifier.%s;',params{i}));
            eval(sprintf('[%s, direction.%s] = getnextpoint(%s,direction.%s,modelpoints.%s,match,sigmas.%s,dists.%s,densityfunction.%s,nextpointweight.%s,parameters);',params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i},params{i}))
        end
    end
end