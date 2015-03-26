function vox = displaysingleflame(flame,outpath)

h = gcf;

if(~iscell(flame))
    t = flame;
    flame = cell(1);
    flame{1} = t;
end

for i = length(flame):-1:1
    if(length(flame) > 1)
        progressbar((length(flame)-i)/length(flame))
    end
    thisflame = flame{i};
%     vox{i} = newrenderskeleton( thisflame, 128, [], true );
%     vox{i} = vox{i}./5;
    vox{i} = renderflame(thisflame,128);
    vox{i}(:) = min(vox{i}(:),1);
    alpha{i} = [];
end
for i = 1:length(vox)
    figure(h);clf;
    vox{i} = gfilter(vox{i},[0.8 0.8 0.8 0]);
    alpha{i} = sum(vox{i},4)./3;
    vol3d('cdata',vox{i},'alpha',alpha{i});
    if(i == 1)
        set(gca,'Color',[0 0 0]);
        set(gcf,'Color',[0 0 0]);
        set(gca,'XtickLabel',[],'YtickLabel',[]);
        set(gca,'visible','off')
        set(gcf, 'InvertHardcopy', 'off');
    end
    view(80,0);
    if(exist('outpath','var'))
        saveas(gcf,sprintf(outpath,i),'png');
    elseif(length(vox) > 1)
        pause
    end
end
if(nargout > 0 && length(vox)==1)
    vox = vox{1};
end

end