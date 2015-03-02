function [ y ] = getS( mode,x,m,n )
%GETS Summary of this function goes here
%   Detailed explanation goes here

% load S.mat
global S

if mode == 0
    y = 0;
elseif mode == 1
    y = S*x;
elseif mode == 2
    y = S'*x;
end

end

