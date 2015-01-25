function y = densityestimate( test, X, H )
%DENSITYESTIMATE Gives a density estimate at points in test from original
%density X given kernel H
%   Detailed explanation goes here

num = size(test,2);
N = size(X,2);
y = zeros(1,num);
for i = 1:num
    summ = 0;
    for j = 1:N
        summ = summ + det(H)^-0.5*mvnpdf((H^-0.5)*(test(:,i)-X(:,j)));
    end  
    summ = summ./N;
    y(i) = summ;
end

end

