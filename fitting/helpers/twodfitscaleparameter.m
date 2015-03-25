function [ s ] = twodfitscaleparameter(xdata, ydata, mu, sig, uest, sest, m)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


nn = 1/1024;
xdata = (round(xdata./nn)).*nn;

ydata = arrayfun(@(x)mean(ydata(xdata==x)),unique(xdata))';

bins = nn:nn:0.2;
ybins = zeros(1,length(bins));
ybins(1:length(ydata)) = ydata;

xdata = bins;
ydata = ybins;

step = @(m,u,x)(heaviside(m-x).*u);

modelfun = @(u,s,x)(rotandsum(s.*(...%TG%
    (step(m,1-u,x)).*(...
        exp(-((x-mu).^2)./(2*sig.^2)))...
    + step(m,u,x))));


beta0 = [uest sest];
lower = [0 0];
upper = [1 0.5];


[f,gof] = fit(xdata',ydata',fittype(modelfun),'Lower',lower,'Upper',upper,'StartPoint',beta0,'Display','off');%,'Robust','Bisquare');

s.s = f.s;
s.u = f.u;

end

