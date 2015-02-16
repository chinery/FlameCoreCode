function [ output ] = rotandsum( input )
%ROTANDSUM numerically rotates and sums a 1d input
%   assume input is x = 0 -> end

trans = false;
if(size(input,1) > 1)
    input = input';
    trans = true;
end

% save('/Users/andy/Dropbox/Work/PhD/FitToImage/buggy.mat','input');

% trim zeros
len = length(input);
trim = find(abs(input) > 0,1,'last');
input(trim+1:len) = [];


[colrange,rowrange] = meshgrid(1:length(input),1:length(input)*2-1);
[colrange,rowrange] = meshgrid(1:length(input),1:length(input));

mesh = [rowrange(:)'; colrange(:)'];

dis = vnorm(bsxfun(@minus,mesh,[length(input);1]));

values = interp1(0:length(input)-1,input,dis,'cubic',0);

grid = reshape(values,length(input),length(input));
grid = [grid ; flipud(grid)];

output = sum(grid,1);

output(trim+1:len) = 0;

output = output./length(input);

if(trans)
    output = output';
end

% figure(10);hold on;plot(linspace(0,0.2,length(output)),output,'-');

end

