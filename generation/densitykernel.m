function [ H ] = densitykernel( X )
%DENSITYKERNEL Summary of this function goes here
%   Detailed explanation goes here

[dim N] = size(X);
H = zeros(dim,dim);
for i = 1:dim
    H(i,i) = (((4/(dim+2))^(1/(dim+4)))*(N^(-1/(dim+4)))*std(X(i,:)))^2;
end

end

