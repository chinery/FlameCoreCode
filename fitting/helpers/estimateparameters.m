function p = estimateparameters(xdata,ydata,prior,tight,col,num)

% figure(1);
% plot(xdata,ydata,'.');

if(nargin < 4)
    tight = false;
end

if(nargin < 3)
    prior = [];
end

nn = 1/1024;
xdata = (round(xdata./nn)).*nn;

ydata = arrayfun(@(x)mean(ydata(xdata==x)),unique(xdata))';

orglength = length(ydata);

global knn

if orglength > 204
    bins = nn:nn:(max(xdata)+0.1);
    ybins = zeros(1,length(bins));
    ybins(1:length(ydata)) = ydata;
    
    label = predict(knn,ybins(1:204));
else
    bins = nn:nn:0.2;
    ybins = zeros(1,length(bins));
    ybins(1:length(ydata)) = ydata;

    label = predict(knn,ybins);
end



mest = bins(length(ydata));
if(label == 1 || label == 4 || label == 2)
    sest = max(ydata);
    
    sigest = 0.1;
    
    if(label == 2)
        [~,ix] = max(ydata);
    else
        [~,ix] = max(ydata.*linspace(1,2,length(ydata)));
    end
    muest = bins(ix);

    lower = [0 0 0 0 0];
    upper = [mest*2 mest 1 0.5 1];
elseif(label == 3)
    sest = median(ydata);
    
    sigest = mest/2;
    
    muest = mest/2;
    
    lower = [-mest/2 mest/4 0 0 0];
    upper = [mest*2 mest 1 0.5 1];
end
uest = sest/min(ydata);


weights = ones(1,length(bins))*0.1;
weights(1:length(ydata)) = 1;

xdata = bins;
ydata = ybins;


% figure(2);
% plot(xdata,ydata,'.');
% pause;

% figure(1)
% clf; hold on;
% plot(xdata,ydata,'.')

step = @(m,u,x)(heaviside(m-x).*u);

            
if(isempty(prior))
%     modelfun = @(l,k,u,s,x)(s.*(... 
%         (1-u).*(...
%             ((x./l).^(k-1)).*exp(-(x./l).^k))...
%         + u));
    
    modelfun = @(mu,sig,u,s,x)(s.*(... %TG%
    (1-u).*(...
        exp(-((x-mu).^2)/(2*sig^2)))...
    + u));
% modelfun = @(mu,sig,u,s,x)(s.*(... %TG%
%     exp(-((x-mu).^2)./(2*sig^2)))...
% + u);



%     beta0 = [0.5 1 ydata(1) 0.1]; 
    beta0 = [muest sigest uest sest]; %TG%
    
    l = lower([1 2 3 5]);
    u = upper([1 2 3 5]);

%     f = fit(xdata',ydata',fittype(modelfun),'Weights',weights,'Lower',[0 0 0 0],'Upper',[1 20 1 1],'StartPoint',beta0,'Display','off');%,'Robust','Bisquare');
    f = fit(xdata',ydata',fittype(modelfun),'Weights',weights,'Lower',l,'Upper',u,'StartPoint',beta0,'Display','off');
    
%     beta0 = [f.l f.k f.u mest f.s];
    beta0 = [f.mu f.sig f.u mest f.s];
else
    beta0 = [prior(1) prior(2) uest mest sest];
%     beta0 = prior;
end

% modelfun = @(l,k,u,m,s,x)(s.*(... 
%     (step(m,1-u,x)).*(...
%         ((x./l).^(k-1)).*exp(-(x./l).^k))...
%     + step(m,u,x)));

modelfun = @(mu,sig,u,m,s,x)(s.*(...%TG%
    (step(m,1-u,x)).*(...
        exp(-((x-mu).^2)./(2*sig.^2)))...
    + step(m,u,x)));
% modelfun = @(mu,sig,u,m,s,x)(step(m,s,x).*(... %TG%
%     exp(-((x-mu).^2)./(2*sig.^2)))...
% + step(m,u,x));

% % % if(~isempty(prior) && tight)
% % % %     lower = max([0 0 0 0 0],prior-abs(0.2*[1 0.5 1 0.5 1]));
% % % %     upper = min([1 0.5 1 0.5 1],prior+abs(0.2*[1 0.5 1 0.5 1]));
% % %     lower = [0 0 0 0 0];
% % %     upper = [1 0.01 1 0.5 1];
% % %     lower(1:3) = max([0 0 0],prior(1:3)-[0.003 0.003 0.05]);
% % %     upper(1:3) = min([1 0.5 1],prior(1:3)+[0.003 0.003 0.05]);
% % % else
% % %     lower = [0 0 0 0 0];
% % %     upper = [mest*2 mest 1 0.5 1];
% % % end

% f = fit(xdata',ydata',fittype(modelfun),'Lower',[0 0 0 0 0],'Upper',[1 20 1 0.5 1],'StartPoint',beta0,'Display','off');%,'Robust','Bisquare');
[f,gof] = fit(xdata',ydata',fittype(modelfun),'Lower',lower,'Upper',upper,'StartPoint',beta0,'Display','off');%,'Robust','Bisquare');

% if(tight)
%     disp(prior(1:2))
% end

% % if(~isempty(prior) && f.sig < 0.5*prior(2) && f.mu < f.m && prior(1) < prior(4))
% %    lower(2) = prior(2); 
% %    [f,gof] = fit(xdata',ydata',fittype(modelfun),'Lower',lower,'Upper',upper,'StartPoint',beta0,'Display','off');%,'Robust','Bisquare'); 
% % end

% % if(f.mu+2*f.sig < f.m || f.sig < 0.015)
% %     beta0 = [0 0.05 beta0(3:5)];
% %     [newf1, newgof] = fit(xdata',ydata',fittype(modelfun),'Lower',lower,'Upper',upper,'StartPoint',beta0,'Display','off');%,'Robust','Bisquare');
% %     if(newgof.rmse < gof.rmse || (f.sig < 0.015 && f.s > 0.15))
% %         f = newf1;
% %         gof = newgof;
% %     end
% %     beta0 = [f.m/2 0.05 beta0(3:5)];
% %     [newf2, newgof] = fit(xdata',ydata',fittype(modelfun),'Lower',lower,'Upper',upper,'StartPoint',beta0,'Display','off');%,'Robust','Bisquare');
% %     if(newgof.rmse < gof.rmse)
% %         f = newf2;
% %     end
% % end
% 
% if(tight)
% figure(1)
% clf; hold on;
% x = linspace(0,0.2,100);
% plot(xdata,ydata,'.')
% axis([0 0.2 0 0.45]);
% if(exist('newf1','var'))
% plot(x,newf1(x),'-g');
% plot(x,newf2(x),'-b');
% end
% plot(x,f(x),'-r');
% % % % export_fig(sprintf('C:/Users/Administrator/Dropbox/Work/PhD/MyData2014/candle1-2/modeltest-l/%i-%01i.png',col,num));
% pause(0.1);
% pause
% end

% i = input('class?');
% 
% global classify
% classify{i}(end+1,:) = ydata;

p.mu = f.mu;
p.sig = f.sig;
p.u = f.u;
p.m = f.m;
p.s = f.s;
p.rmse = gof.rmse;