function [ p ] = advanceddensityestimate( test, pdf )
%ADVANCEDDENSITYESTIMATE Summary of this function goes here
%   Detailed explanation goes here

p = evaluatePointsUnderPdf(pdf, test);

end

