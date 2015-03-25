% to work out how the licks move from flame A to B
function match = matchflames( frompoints, fromlabel, topoints, tolabel )
%SIMPLEMATCH Summary of this function goes here
%   Detailed explanation goes here


points1 = frompoints;
points2 = topoints;

label1 = fromlabel;
label2 = tolabel;

num1 = max(label1)-1;
num2 = max(label2)-1;

prob = zeros(num1,num2);
for i = 1:num1
    part = points1(:,label1==i);
    if (length(part) < 10)
        continue;
    end
    mu = mean(part,2);
    sig = cov(part');
    
    mvn1 = mvn_new(sig,mu);
    
    for j = 1:num2
        part = points2(:,label2==j);
        if (length(part) < 10)
            continue;
        end
        mu2 = mean(part,2);
        sig2 = cov(part');

        mvn2 = mvn_new(sig2,mu2);     
        
%         % weight based on angle up-ness
%         v = mu2 - mu;
%         angle = atan2(norm(cross(v,[0 0 1])),dot(v,[0 0 1]));
%         prior = 1;%-(angle/pi);
        
        dist = mvn_div_kl(mvn1,mvn2)/20;
%         dist = norm(mu-mu2)/30;
%         probs(i,j) = prior*exp(-1*dist);
        prob(i,j) = exp(-1*dist);
    end
    
end

prob(:,end+1) = 0.1;
prob = bsxfun(@rdivide,prob,sum(prob,2));

if(true)%num1 == num2 || all(max(prob,[],2) > 0.5))
    [val, mat] = max(prob,[],1);
else
    [val, mat] = max(prob,[],1);
    figure(9);clf;hold on;
    for j = 1:num1
        subplot(num1,2,sub2ind([2 num1],1,j)); hold on;
        for k = 1:num1
            part = points1(:,label1==k);
            
            if(mat(j) == k)
                plot3(part(1,:),part(2,:),part(3,:),'.g');
            else
                plot3(part(1,:),part(2,:),part(3,:),'.b');
            end
        end
        view(20,15);
        
        subplot(num1,2,sub2ind([2 num1],2,j)); hold on;
        for k = 1:num2
            part = points2(:,label2==k);
            
            if(j == k)
                plot3(part(1,:),part(2,:),part(3,:),'.g');
            else
                plot3(part(1,:),part(2,:),part(3,:),'.b');
            end
        end
        view(20,15);
    end
    def = { num2str(1) };
    options.WindowStyle = 'normal';
    x = inputdlg('This okay? \n 1: Yes 0: No', 'Uncertain matching', 1, def,options);
    if(isempty(x) || ~str2double(x{:}))
        error('abort!');
    end
end


tomask = true(1,max(tolabel));
tomask(end) = false;
frommask = true(1,max(fromlabel));
frommask(end) = false;

figure(9); clf; hold on;
toooo = find(tomask);
for jix = 1:sum(tomask)
    j = toooo(jix);
    subplot(sum(tomask),2,sub2ind([2 sum(tomask)],1,jix)); hold on;
    for k = 1:num1
        part = points1(:,label1==k);
        
        if(mat(j) == k)
            plot3(part(1,:),part(2,:),part(3,:),'.g','markersize',1);
        elseif(~frommask(k))
%             plot3(part(1,:),part(2,:),part(3,:),'.r','markersize',1);
        else
            plot3(part(1,:),part(2,:),part(3,:),'.','color',[0.7 0.7 0.7],'markersize',1);
        end
    end
    axis equal
    view(20,15);
    
    subplot(sum(tomask),2,sub2ind([2 sum(tomask)],2,jix));  hold on;
    for k = 1:num2
        part = points2(:,label2==k);
        
        if(j == k)
            plot3(part(1,:),part(2,:),part(3,:),'.g','markersize',1);
        elseif(~tomask(k))
%             plot3(part(1,:),part(2,:),part(3,:),'.r','markersize',1);
        else
            plot3(part(1,:),part(2,:),part(3,:),'.','color',[0.7 0.7 0.7],'markersize',1);
        end
    end
    axis equal
    view(20,15);
end


left = mat(1:num2);
right = 1:num2;

match = [left; right];

% if(any(~frommask(match(1,:))))
%     tomask(~frommask(match(1,:))) = [];
%     val(~frommask(match(1,:))) = [];
%     match(:,~frommask(match(1,:))) = [];
% end

if(any(val == 0))
    match(:,val == 0) = [];
end

end

