addpath([pwd '/common/']);
addpath([pwd '/external/']);

folder = [pwd '/external/'];
sub = dir(folder);

for i = 3:length(sub)
    if(sub(i).isdir)
       addpath([folder sub(i).name '/']);
    end
end